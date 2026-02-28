/*
 * Cessna++ AWG (v35) – Phase-Segmented, Servo-Tracking Arbitrary Waveform Generator
 * --------------------------------------------------------------------------------
 * Design decisions implemented here (as per conversation):
 * - One-cycle segmentation into 16 bins.
 * - Bin value = mean (average) of input samples within each bin.
 * - Bins computed once per completed cycle.
 * - Servo convergence for bin updates: A[i] += alpha*(new[i]-A[i]).
 * - Per-bin deviation from Sweet Sixteen via I2C mappings: S[i] = A[i] + D*X[i].
 *   (Expect user to enable "Symmetric" mapping on the NT for bipolar control.)
 * - Deviations applied pre-interpolation.
 * - Interpolation morph: Linear <-> Curved (shape-preserving-ish cubic Hermite).
 * - Cycle detection: rising zero-crossing with auto hysteresis threshold.
 * - Freeze: stops bin updates (alpha forced to 0), but cycle detection keeps running.
 * - Capture: only when LOCKED, capture next full cycle and then Freeze ON.
 * - Outputs: OUT1 = resynth audio, OUT2 = dry input.
 * - Display: overlays one-cycle input and resynth waveforms; bottom ticks 1–16 at bin centers;
 *   tick brightens briefly when corresponding X[i] changes.
 */

#include <math.h>
#include <string.h>
#include <new>
#include <distingnt/api.h>

#ifndef M_PI_F
#define M_PI_F 3.14159265358979323846f
#endif

extern "C" int* __errno(void) {
    static int errno_val = 0;
    return &errno_val;
}

// ----------------------------
// Constants / configuration
// ----------------------------
static constexpr int   kBins              = 16;
static constexpr int   kDisplayPoints     = 128;   // fixed points per cycle for display (phase-normalised)
static constexpr float kDefaultSampleRate = 48000.0f;

// Cycle detection guard rails (Hz)
static constexpr float kMinHz = 20.0f;
static constexpr float kMaxHz = 5000.0f; // conservative; can be raised later

// Auto-threshold parameters (revised in processCycleDetection; kept here for reference)
static constexpr float kEnvBeta        = 0.01f;    // envelope smoothing (per sample) — was 0.0015f, increased for faster threshold adaptation

// Max samples to keep for one captured input cycle
// (at 48kHz and 20Hz minimum, one cycle is 2400 samples; add a little margin)
static constexpr int kMaxCycleSamps = (int)(kDefaultSampleRate / kMinHz) + 64;

// Minimum cycle length before bin averaging is considered reliable.
// Below this (approx 1kHz at 48kHz), bin updates are skipped to avoid
// corrupting A[] with noisy single-sample averages on high notes.
static constexpr int kMinFinalizeSamps = 48;

// ----------------------------
// Helpers
// ----------------------------
static inline float clampf(float x, float lo, float hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

static inline float lerpf(float a, float b, float t) {
    return a + (b - a) * t;
}

// ---- Forward declarations (MUST appear before construct/step) ----
static inline float evalLinear(const float *bins, float phase01);
static bool processCycleDetection(struct _CessnaAwg_DTC *d, float x);
static void finalizeCycle(struct _CessnaAwg_DTC *d, int cycleLen);
static void seqBuildPattern(struct _CessnaAwg_DTC *d);
static void computeModX(struct _CessnaAwg_DTC *d);

// Simple DC blocker (one-pole HPF)
struct DcBlocker {
    float x1 = 0.0f;
    float y1 = 0.0f;
    float a  = 0.995f; // closer to 1 = lower cutoff

    inline float process(float x) {
        // y[n] = x[n] - x[n-1] + a*y[n-1]
        float y = x - x1 + a * y1;
        x1 = x;
        y1 = y;
        return y;
    }
};

// Compute shape-preserving-ish cubic Hermite for periodic control points.
struct PeriodicHermite16 {
    float m[kBins] = {0}; // slopes at points

    void computeSlopes(const float *p) {
        float d[kBins];
        for (int i=0;i<kBins;i++) {
            float p0 = p[i];
            float p1 = p[(i+1)&15];
            d[i] = (p1 - p0);
        }

        for (int i=0;i<kBins;i++) {
            float dm = d[(i-1)&15];
            float dp = d[i];
            m[i] = 0.5f * (dm + dp);
        }

        // Pass 1: zero slopes where sign conflicts with the chord (Fritsch-Carlson).
        // Must be a separate pass so we don't corrupt a bin's slope before its
        // neighbour reads it as input.
        for (int i=0;i<kBins;i++) {
            float dp = d[i];
            if (dp == 0.0f) {
                m[i] = 0.0f;
                m[(i+1)&15] = 0.0f;
                continue;
            }
            if ((m[i] / dp) < 0.0f)        m[i]        = 0.0f;
            if ((m[(i+1)&15] / dp) < 0.0f) m[(i+1)&15] = 0.0f;
        }

        // Pass 2: norm-clamp |m| <= 3*|dp| for shape preservation.
        for (int i=0;i<kBins;i++) {
            float dp = d[i];
            if (dp == 0.0f) continue;
            float r1 = m[i] / dp;
            float r2 = m[(i+1)&15] / dp;
            float sumsq = r1*r1 + r2*r2;
            if (sumsq > 9.0f) {
                float t = 3.0f / sqrtf(sumsq);
                m[i]        *= t;
                m[(i+1)&15] *= t;
            }
        }
    }

    inline float eval(const float *p, float phase01) const {
        float x = phase01 * (float)kBins;
        int i = (int)floorf(x);
        float t = x - (float)i;
        i &= 15;
        int i1 = (i+1) & 15;

        float p0 = p[i];
        float p1 = p[i1];
        float m0 = m[i];
        float m1 = m[i1];

        float t2 = t*t;
        float t3 = t2*t;
        float h00 =  2.0f*t3 - 3.0f*t2 + 1.0f;
        float h10 =        t3 - 2.0f*t2 + t;
        float h01 = -2.0f*t3 + 3.0f*t2;
        float h11 =        t3 -       t2;

        return h00*p0 + h10*m0 + h01*p1 + h11*m1;
    }
};

// ----------------------------
// DTC (fast memory)
// ----------------------------
struct _CessnaAwg_DTC {
    float sampleRate;

    float depthD;
    float adaptAlpha;
    float curve;
    float lockStrength;
    float trackStrength;    // 0=free-run, 1=full frequency tracking
    bool  freeze;
    bool  pendingCapture;

    // Sync mode: 0=Edge, 1=N-before-M, 2=PLL
    int   syncMode;

    // N-before-M state
    int   nbmN;             // lead/lag register length (derived from Lock%)
    int   nbmM;             // total input sum register length (= 2*N, as per Cessna)
    int   nbmLeadCount;
    int   nbmLagCount;
    int   nbmTotalCount;

    // PLL state
    float pllPhaseErr;      // smoothed phase error (radians equivalent, 0..1)
    float pllIntegrator;    // loop filter integrator
    float pllBandwidth;     // loop bandwidth coefficient (derived from Lock%)

    float X[kBins];

    DcBlocker dc;
    float env;
    float thresh;
    bool  wasBelowNeg;
    int   samplesSinceEdge;
    int   lastPeriodSamples;
    float smoothPeriod;
    bool  locked;
    bool  hardJump;         // true when pitch jump detected — skip finalizeCycle

    int   prevPeriodSamples;
    int   stableCount;

    float cycleBuf[kMaxCycleSamps];
    int   cycleWrite;
    int   minPeriodSamps;
    int   maxPeriodSamps;

    float A[kBins];
    float newA[kBins];
    float S[kBins];

    PeriodicHermite16 herm;

    float phase;
    float phaseInc;
    float phaseIncTarget;   // slew target for phaseInc
    float slewCoeff;        // per-sample slew coefficient (set in construct)

    // --- Internal modulation state ---
    int   modMode;          // 0=Off 1=Ripple 2=Basis 3=Sequencer
    float modDepth;         // 0..1
    float modRate;          // Hz
    float modSpread;        // 0..1 (Ripple)
    float lfoPhase;         // 0..1, advances each sample
    float lfoPhaseInc;      // per-sample increment

    // Basis mode
    float basisStatic[4];   // A,B,C,D static coefficients -1..1
    float basisMix;         // 0=static, 1=LFO

    // Sequencer mode
    int   seqSpeed;         // cycles per step
    float seqSmooth;        // 0..1 slew
    float seqDrift;         // 0..1 random walk per cycle
    int   seqPattern;       // 0=Ramp 1=Triangle 2=Random 3=Alternating
    float seqTarget[kBins]; // target values (pattern)
    float seqCurrent[kBins];// slewed current values
    int   seqCycleCount;    // cycles since last step
    uint32_t seqRng;        // simple LCG state for random

    // Combined modulation output (computed each sample in Ripple/Basis, each cycle in Seq)
    float modX[kBins];      // internal mod contribution, -1..1

    int   lastChangedParam; // most recently edited parameter, for display
    bool  needsGrayOutUpdate; // deferred grey-out on first step()

    float dispIn[kDisplayPoints];
    float dispOut[kDisplayPoints];
    bool  dispInit;

    _CessnaAwg_DTC() {
        sampleRate = kDefaultSampleRate;

        depthD = 0.0f;
        adaptAlpha = 1.0f;
        curve = 0.0f;
        lockStrength = 1.0f;
        trackStrength = 1.0f;
        freeze = false;
        pendingCapture = false;

        syncMode = 0;

        nbmN = 4;
        nbmM = 8;
        nbmLeadCount  = 0;
        nbmLagCount   = 0;
        nbmTotalCount = 0;

        pllPhaseErr   = 0.0f;
        pllIntegrator = 0.0f;
        pllBandwidth  = 0.05f;

        env = 0.0f;
        thresh = 0.0f;
        wasBelowNeg = false;
        samplesSinceEdge = 0;
        lastPeriodSamples = (int)(kDefaultSampleRate / 110.0f);
        smoothPeriod = (float)lastPeriodSamples;
        locked = false;
        hardJump = false;

        prevPeriodSamples = lastPeriodSamples;
        stableCount = 0;

        cycleWrite = 0;
        minPeriodSamps = (int)floorf(sampleRate / kMaxHz);
        maxPeriodSamps = (int)ceilf(sampleRate / kMinHz);
        if (minPeriodSamps < 4)  minPeriodSamps = 4;
        if (minPeriodSamps > 32) minPeriodSamps = 32;  // hard cap
        if (maxPeriodSamps > kMaxCycleSamps) maxPeriodSamps = kMaxCycleSamps;

        phase = 0.0f;
        phaseInc = 1.0f / smoothPeriod;
        phaseIncTarget = phaseInc;
        slewCoeff = 0.0f; // set properly in construct()

        modMode     = 0;
        modDepth    = 0.5f;
        modRate     = 0.1f;
        modSpread   = 1.0f;
        lfoPhase    = 0.0f;
        lfoPhaseInc = modRate / kDefaultSampleRate;

        for (int i=0;i<4;i++) basisStatic[i] = (i==0) ? 1.0f : 0.0f;
        basisMix = 0.0f;

        seqSpeed      = 4;
        seqSmooth     = 0.5f;
        seqDrift      = 0.2f;
        seqPattern    = 0;
        seqCycleCount = 0;
        seqRng        = 0x12345678u;
        for (int i=0;i<kBins;i++) {
            seqTarget[i]  = 0.0f;
            seqCurrent[i] = 0.0f;
            modX[i]       = 0.0f;
        }

        lastChangedParam = 0;
        needsGrayOutUpdate = true;

        dispInit = false;

        for (int i=0;i<kBins;i++) {
            X[i] = 0.0f;

            float ph = ((float)i + 0.5f) / (float)kBins;
            A[i] = sinf(2.0f*M_PI_F*ph) * 5.0f;
            newA[i] = A[i];
            S[i] = A[i];
        }

        for (int i=0;i<kMaxCycleSamps;i++) cycleBuf[i] = 0.0f;
        for (int i=0;i<kDisplayPoints;i++) { dispIn[i]=0.0f; dispOut[i]=0.0f; }
    }
};

// ----------------------------
// Algorithm object
// ----------------------------
struct _CessnaAwg : public _NT_algorithm {
    _CessnaAwg_DTC *dtc;
    _CessnaAwg(_CessnaAwg_DTC *d) : dtc(d) {}
};

// ----------------------------
// Parameters
// ----------------------------
enum {
    kParamInput = 0,
    kParamOut1, kParamOut1Mode,
    kParamOut2, kParamOut2Mode,
    kParamOut3, kParamOut3Mode,

    kParamDepth,
    kParamAdapt,
    kParamCurve,
    kParamSyncMode,     // 0=Edge 1=N-before-M 2=PLL
    kParamLock,
    kParamTrack,
    kParamPllBw,        // PLL bandwidth 1..100 (temp tuning parameter)
    kParamFreeze,
    kParamCapture,

    kParamX1,
    kParamX16 = kParamX1 + 15,

    // Mod page
    kParamModMode,      // 0=Off 1=Ripple 2=Basis 3=Sequencer
    kParamModRate,      // 1..1000 → 0.01..10 Hz (stored as integer centimHz)
    kParamModDepth,     // 0..100%
    kParamModSpread,    // 0..100% (Ripple)
    kParamBasisA,       // -100..100% (Basis static)
    kParamBasisB,
    kParamBasisC,
    kParamBasisD,
    kParamBasisMix,     // 0..100% (static→LFO blend)
    kParamSeqSpeed,     // 1..32 cycles per step (Sequencer)
    kParamSeqSmooth,    // 0..100% slew (Sequencer)
    kParamSeqDrift,     // 0..100% random walk (Sequencer)
    kParamSeqPattern,   // 0=Ramp 1=Triangle 2=Random 3=Alternating (Sequencer)
};

static const char* onOffStrings[]    = {"Off", "On", nullptr};
static const char* syncModeStrings[] = {"Edge", "N-before-M", "PLL", nullptr};
static const char* modModeStrings[]  = {"Off", "Ripple", "Basis", "Sequencer", nullptr};
static const char* seqPatStrings[]   = {"Ramp", "Triangle", "Random", "Alternating", nullptr};
static constexpr int kNumParams = (kParamSeqPattern + 1);

// ----------------------------
// Parameters (new API)
// ----------------------------
static const _NT_parameter g_parameters[kNumParams] = {
    NT_PARAMETER_AUDIO_INPUT("in", 0, 1)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("resynth", 0, 13)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("dry", 0, 14)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("error", 0, 0)

    { .name="Depth",     .min=0,   .max=100, .def=0,   .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Adapt",     .min=0,   .max=100, .def=100, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Curve",     .min=0,   .max=100, .def=0,   .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Sync Mode", .min=0,   .max=2,   .def=0,   .unit=kNT_unitEnum,    .scaling=kNT_scalingNone, .enumStrings=syncModeStrings },
    { .name="Lock",      .min=0,   .max=100, .def=100, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Track",     .min=0,   .max=100, .def=100, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="PLL BW",    .min=1,   .max=100, .def=5,   .unit=kNT_unitNone,    .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="freeze",    .min=0,   .max=1,   .def=0,   .unit=kNT_unitEnum,    .scaling=kNT_scalingNone, .enumStrings=onOffStrings },
    { .name="capture",   .min=0,   .max=1,   .def=0,   .unit=kNT_unitEnum,    .scaling=kNT_scalingNone, .enumStrings=onOffStrings },

    { .name="Offset 1",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 2",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 3",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 4",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 5",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 6",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 7",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 8",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 9",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 10", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 11", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 12", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 13", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 14", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 15", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Offset 16", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },

    // Mod page
    { .name="Mod Mode",  .min=0,.max=3,   .def=0,  .unit=kNT_unitEnum,    .scaling=kNT_scalingNone, .enumStrings=modModeStrings },
    { .name="Rate",      .min=1,.max=1000,.def=10, .unit=kNT_unitNone,    .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Mod Depth", .min=0,.max=100, .def=50, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Spread",    .min=0,.max=100, .def=100,.unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Basis A",   .min=-100,.max=100,.def=100,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Basis B",   .min=-100,.max=100,.def=0,  .unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Basis C",   .min=-100,.max=100,.def=0,  .unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Basis D",   .min=-100,.max=100,.def=0,  .unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Basis Mix", .min=0,.max=100,  .def=0,  .unit=kNT_unitPercent, .scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Seq Speed", .min=1,.max=32,   .def=4,  .unit=kNT_unitNone,    .scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Seq Smooth",.min=0,.max=100,  .def=50, .unit=kNT_unitPercent, .scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Seq Drift", .min=0,.max=100,  .def=20, .unit=kNT_unitPercent, .scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Seq Pattern",.min=0,.max=3,   .def=0,  .unit=kNT_unitEnum,    .scaling=kNT_scalingNone,.enumStrings=seqPatStrings },
};

static const uint8_t g_pageRoutingIdx[] = {
    (uint8_t)kParamInput,
    (uint8_t)kParamOut1, (uint8_t)kParamOut1Mode,
    (uint8_t)kParamOut2, (uint8_t)kParamOut2Mode,
    (uint8_t)kParamOut3, (uint8_t)kParamOut3Mode,
};

static const uint8_t g_pageMainIdx[] = {
    (uint8_t)kParamDepth,
    (uint8_t)kParamAdapt,
    (uint8_t)kParamCurve,
    (uint8_t)kParamSyncMode,
    (uint8_t)kParamLock,
    (uint8_t)kParamTrack,
    (uint8_t)kParamPllBw,
    (uint8_t)kParamFreeze,
    (uint8_t)kParamCapture,
};

static const uint8_t g_pageXIdx[] = {
    (uint8_t)kParamX1,  (uint8_t)(kParamX1+1),  (uint8_t)(kParamX1+2),  (uint8_t)(kParamX1+3),
    (uint8_t)(kParamX1+4),  (uint8_t)(kParamX1+5),  (uint8_t)(kParamX1+6),  (uint8_t)(kParamX1+7),
    (uint8_t)(kParamX1+8),  (uint8_t)(kParamX1+9),  (uint8_t)(kParamX1+10), (uint8_t)(kParamX1+11),
    (uint8_t)(kParamX1+12), (uint8_t)(kParamX1+13), (uint8_t)(kParamX1+14), (uint8_t)(kParamX1+15),
};

static const uint8_t g_pageModIdx[] = {
    (uint8_t)kParamModMode,
    (uint8_t)kParamModRate,
    (uint8_t)kParamModDepth,
    (uint8_t)kParamModSpread,
    (uint8_t)kParamBasisA,
    (uint8_t)kParamBasisB,
    (uint8_t)kParamBasisC,
    (uint8_t)kParamBasisD,
    (uint8_t)kParamBasisMix,
    (uint8_t)kParamSeqSpeed,
    (uint8_t)kParamSeqSmooth,
    (uint8_t)kParamSeqDrift,
    (uint8_t)kParamSeqPattern,
};

static const _NT_parameterPage g_pages_arr[] = {
    { .name="Routing", .numParams=(uint8_t)ARRAY_SIZE(g_pageRoutingIdx), .group=0, .unused={0,0}, .params=g_pageRoutingIdx },
    { .name="Main",    .numParams=(uint8_t)ARRAY_SIZE(g_pageMainIdx),    .group=0, .unused={0,0}, .params=g_pageMainIdx },
    { .name="Offsets", .numParams=(uint8_t)ARRAY_SIZE(g_pageXIdx),       .group=1, .unused={0,0}, .params=g_pageXIdx },
    { .name="Mod",     .numParams=(uint8_t)ARRAY_SIZE(g_pageModIdx),     .group=0, .unused={0,0}, .params=g_pageModIdx },
};

static const _NT_parameterPages g_pages = {
    .numPages = ARRAY_SIZE(g_pages_arr),
    .pages = g_pages_arr,
};

// ----------------------------
// Requirements / construct (new API)
// ----------------------------
static void calculateRequirements(_NT_algorithmRequirements &req, const int32_t* specifications) {
    (void)specifications;
    req.numParameters = kNumParams;
    req.sram = sizeof(_CessnaAwg);
    req.dram = 0;
    req.dtc  = sizeof(_CessnaAwg_DTC);
    req.itc  = 0;
}

static _NT_algorithm* construct(const _NT_algorithmMemoryPtrs& ptrs,
                               const _NT_algorithmRequirements& req,
                               const int32_t* specifications) {
    (void)req;
    (void)specifications;

    auto *d = new(ptrs.dtc) _CessnaAwg_DTC();
    auto *self = new(ptrs.sram) _CessnaAwg(d);

    self->parameters = g_parameters;
    self->parameterPages = &g_pages;

    // Guard sample rate — NT_globals.sampleRate may be uninitialised or out of range.
    // Accept only plausible audio rates (8kHz–192kHz), otherwise fall back to default.
    float sr = (float)NT_globals.sampleRate;
    if (sr < 8000.0f || sr > 192000.0f) sr = kDefaultSampleRate;
    d->sampleRate = sr;

    // Hard-code sane period bounds rather than deriving from sample rate,
    // since NT_globals.sampleRate may not be trustworthy at construct time.
    // minPeriodSamps: 5000Hz max → 48000/5000 = ~10. Use 8 as a safe floor.
    // maxPeriodSamps: 20Hz min → 48000/20 = 2400.
    d->minPeriodSamps = (int)floorf(sr / kMaxHz);
    d->maxPeriodSamps = (int)ceilf(sr / kMinHz);
    if (d->minPeriodSamps < 4)    d->minPeriodSamps = 4;
    if (d->minPeriodSamps > 32)   d->minPeriodSamps = 32;   // hard cap — never more than 32
    if (d->maxPeriodSamps < 256)  d->maxPeriodSamps = 256;
    if (d->maxPeriodSamps > kMaxCycleSamps) d->maxPeriodSamps = kMaxCycleSamps;

    d->smoothPeriod = (float)d->lastPeriodSamples;
    d->phaseInc = 1.0f / d->smoothPeriod;
    d->phaseIncTarget = d->phaseInc;
    // Slew: ~5ms time constant. coeff = 1 - exp(-1/(sr*0.005))
    // Approximation: coeff = 1/(sr*0.005) for small values
    d->slewCoeff = 1.0f / (sr * 0.005f);
    if (d->slewCoeff > 1.0f) d->slewCoeff = 1.0f;
    d->lfoPhaseInc = d->modRate / d->sampleRate;

    // Seed sequencer with ramp pattern
    seqBuildPattern(d);

    for (int i=0;i<kBins;i++) d->S[i] = d->A[i];
    d->herm.computeSlopes(d->S);
    for (int k=0;k<kDisplayPoints;k++) {
        float ph = (float)k / (float)(kDisplayPoints - 1);
        float yLin = evalLinear(d->S, ph);
        float yCur = d->herm.eval(d->S, ph);
        d->dispOut[k] = lerpf(yLin, yCur, d->curve);
    }
    d->dispInit = true;
    d->needsGrayOutUpdate = true;

    return self;
}

// ----------------------------
// Parameter changes (new API)
// ----------------------------
// ----------------------------
// Grey out irrelevant mod parameters based on current mode
// ----------------------------
static void updateModGrayOut(_NT_algorithm* base) {
    auto *self = (_CessnaAwg*)base;
    if (!self || !self->v) return;

    int mode = self->v[kParamModMode];
    int syncMode = self->v[kParamSyncMode];
    uint32_t idx = (uint32_t)NT_algorithmIndex(self);

    bool isOff  = (mode == 0);
    bool ripple = (mode == 1);
    bool basis  = (mode == 2);
    bool seq    = (mode == 3);

    // PLL BW only relevant in PLL sync mode
    NT_setParameterGrayedOut(idx, kParamPllBw, syncMode != 2);

    // Shared params (Rate, Mod Depth): grey only when Off
    NT_setParameterGrayedOut(idx, kParamModRate,    isOff);
    NT_setParameterGrayedOut(idx, kParamModDepth,   isOff);

    // Ripple-only
    NT_setParameterGrayedOut(idx, kParamModSpread,  !ripple);

    // Basis-only
    NT_setParameterGrayedOut(idx, kParamBasisA,     !basis);
    NT_setParameterGrayedOut(idx, kParamBasisB,     !basis);
    NT_setParameterGrayedOut(idx, kParamBasisC,     !basis);
    NT_setParameterGrayedOut(idx, kParamBasisD,     !basis);
    NT_setParameterGrayedOut(idx, kParamBasisMix,   !basis);

    // Sequencer-only
    NT_setParameterGrayedOut(idx, kParamSeqSpeed,   !seq);
    NT_setParameterGrayedOut(idx, kParamSeqSmooth,  !seq);
    NT_setParameterGrayedOut(idx, kParamSeqDrift,   !seq);
    NT_setParameterGrayedOut(idx, kParamSeqPattern, !seq);
}

static void parameterChanged(_NT_algorithm *base, int p) {
    auto *self = (_CessnaAwg*)base;
    if (!self || !self->dtc) return;
    auto *d = self->dtc;
    if (!self->v) return;

    auto getPct01 = [&](int idx)->float {
        return clampf((float)self->v[idx] / 100.0f, 0.0f, 1.0f);
    };
    auto getPctSigned = [&](int idx)->float {
        return clampf((float)self->v[idx] / 100.0f, -1.0f, 1.0f);
    };

    d->lastChangedParam = p;

    switch (p) {
        case kParamDepth:   d->depthD       = getPct01(kParamDepth); break;
        case kParamAdapt:   d->adaptAlpha   = getPct01(kParamAdapt); break;
        case kParamCurve:   d->curve        = getPct01(kParamCurve); break;
        case kParamSyncMode: {
            d->syncMode = self->v[kParamSyncMode];
            // Reset sync state on mode change
            d->nbmLeadCount  = 0;
            d->nbmLagCount   = 0;
            d->nbmTotalCount = 0;
            d->pllPhaseErr   = 0.0f;
            d->pllIntegrator = 0.0f;
            break;
        }
        case kParamLock: {
            d->lockStrength = getPct01(kParamLock);
            // Derive N-before-M parameters from lock strength:
            // Lock 100% → N=8, M=16 (firm); Lock 0% → N=2, M=4 (loose)
            float ls = d->lockStrength;
            d->nbmN = 2 + (int)(ls * 6.0f);   // 2..8
            d->nbmM = d->nbmN * 2;             // always M=2N (Cessna's recommendation)
            break;
        }
        case kParamTrack:   d->trackStrength = getPct01(kParamTrack); break;
        case kParamPllBw: {
            // param 1..100 → bandwidth 0.001..0.1 (log-ish feel)
            float bw = (float)self->v[kParamPllBw];
            d->pllBandwidth = 0.001f + (bw / 100.0f) * 0.099f;
            break;
        }
        case kParamFreeze:  d->freeze       = (self->v[kParamFreeze] != 0); break;
        case kParamCapture: d->pendingCapture = (self->v[kParamCapture] != 0); break;

        case kParamModMode:  d->modMode  = self->v[kParamModMode]; break;
        case kParamModRate: {
            // param 1..1000 maps to 0.01..10 Hz linearly
            d->modRate     = (float)self->v[kParamModRate] * 0.01f;
            d->lfoPhaseInc = d->modRate / d->sampleRate;
            break;
        }
        case kParamModDepth:  d->modDepth  = getPct01(kParamModDepth); break;
        case kParamModSpread: d->modSpread = getPct01(kParamModSpread); break;
        case kParamBasisA:    d->basisStatic[0] = getPctSigned(kParamBasisA); break;
        case kParamBasisB:    d->basisStatic[1] = getPctSigned(kParamBasisB); break;
        case kParamBasisC:    d->basisStatic[2] = getPctSigned(kParamBasisC); break;
        case kParamBasisD:    d->basisStatic[3] = getPctSigned(kParamBasisD); break;
        case kParamBasisMix:  d->basisMix        = getPct01(kParamBasisMix); break;
        case kParamSeqSpeed:  d->seqSpeed   = self->v[kParamSeqSpeed]; break;
        case kParamSeqSmooth: d->seqSmooth  = getPct01(kParamSeqSmooth); break;
        case kParamSeqDrift:  d->seqDrift   = getPct01(kParamSeqDrift); break;
        case kParamSeqPattern:d->seqPattern = self->v[kParamSeqPattern]; break;
        default:
            if (p >= kParamX1 && p <= kParamX16) {
                int i = p - kParamX1;
                d->X[i] = getPctSigned(p);
            }
            break;
    }
    updateModGrayOut(base);
}

static inline float* busPtr(float* busFrames, int numFrames, int busParamValue) {
    if (!busFrames) return nullptr;
    if (busParamValue <= 0) return nullptr;
    int busIndex = busParamValue - 1;
    if (busIndex < 0 || busIndex >= 28) return nullptr;
    return busFrames + (busIndex * numFrames);
}

static inline float evalLinear(const float *bins, float phase) {
    phase -= floorf(phase);
    float x = phase * (float)kBins;
    int i = (int)floorf(x);
    float f = x - (float)i;
    int i0 = i & (kBins - 1);
    int i1 = (i + 1) & (kBins - 1);
    return lerpf(bins[i0], bins[i1], f);
}

/*
 * Cycle detector — Schmidt trigger with hard hold-off.
 *
 * - Envelope pre-warmed instantly on first samples (no slow ramp from zero).
 * - Threshold raised to 10% of envelope, floored at 0.01 (normalised units).
 * - Hard hold-off of minPeriodSamps after each firing prevents double-triggers
 *   on any waveform shape (square edges, sawtooth reset, etc.).
 * - samplesSinceEdge only resets on an actual firing, not on rejected crossings.
 * - minPeriodSamps check folded into the hold-off guard, so the period range
 *   check only needs an upper bound inside the block.
 */
static bool processCycleDetection(struct _CessnaAwg_DTC *d, float x) {
    // DC-block for envelope estimate only
    float yEnv = d->dc.process(x);
    float ay = fabsf(yEnv);

    // Pre-warm: snap envelope up immediately if signal is louder than current estimate.
    if (ay > d->env) d->env = ay;
    d->env = (1.0f - kEnvBeta) * d->env + kEnvBeta * ay;

    // Adaptive threshold: 10% of RMS envelope, hard floor at 0.01 (normalised).
    const float thresh = fmaxf(0.01f, 0.10f * d->env);
    d->thresh = thresh;

    d->samplesSinceEdge++;

    // Arm when signal goes clearly negative.
    if (x < -thresh) d->wasBelowNeg = true;

    // Fire when: armed, signal clearly positive, AND hold-off elapsed.
    // The hold-off (>= minPeriodSamps) replaces the lower-bound period check,
    // preventing double-fires on fast edges without needing a separate gate.
    if (d->wasBelowNeg && (x > thresh) && (d->samplesSinceEdge >= d->minPeriodSamps)) {
        int period = d->samplesSinceEdge;
        d->samplesSinceEdge = 0;
        d->wasBelowNeg = false;

        if (period <= d->maxPeriodSamps) {
            float rel = fabsf((float)period - (float)d->prevPeriodSamples)
                        / (float)d->prevPeriodSamples;

            if (rel > 0.10f) {
                // Hard pitch jump — update rate but keep phase continuous to avoid click
                d->smoothPeriod = (float)period;
                d->stableCount  = 0;
                d->locked       = false;
                d->hardJump     = true;
            } else {
                // Normal tracking — smooth
                const float beta = 0.10f;
                d->smoothPeriod = (1.0f - beta) * d->smoothPeriod + beta * (float)period;
                if (d->smoothPeriod < 4.0f) d->smoothPeriod = 4.0f;
                if (rel < 0.05f) d->stableCount++; else d->stableCount = 0;
                d->locked   = (d->stableCount >= 3);
                d->hardJump = false;
            }
            d->prevPeriodSamples = period;

            return true;
        } else {
            // period out of range, discard
        }
    }
    return false;
}

// Combine A[i] + scaled offset, keeping S[i] within ±5V.
// The offset pushes within the remaining headroom, so at max offset
// S[i] reaches the rail but never exceeds it.
static __attribute__((noinline)) void combineS(struct _CessnaAwg_DTC *d) {
    float D = d->depthD * 5.0f;
    for (int i = 0; i < kBins; ++i) {
        float offset = D * clampf(d->X[i] + d->modX[i], -1.0f, 1.0f);
        float headroom = (offset >= 0.0f) ? (5.0f - d->A[i]) : (5.0f + d->A[i]);
        headroom = clampf(headroom / 5.0f, 0.0f, 1.0f);
        d->S[i] = d->A[i] + offset * headroom;
    }
}

static void finalizeCycle(struct _CessnaAwg_DTC *d, int cycleLen) {
    if (cycleLen < 8) {
        d->cycleWrite = 0;
        return;
    }

    // Skip bin updates on very short cycles (high notes) — too few samples
    // per bin to get reliable averages. Oscillator keeps running on existing A[].
    if (cycleLen < kMinFinalizeSamps) {
        d->cycleWrite = 0;
        return;
    }

    float sum[kBins] = {0};
    int   cnt[kBins] = {0};

    for (int s = 0; s < cycleLen; ++s) {
        float ph = (float)s / (float)cycleLen;
        int bin = (int)floorf(ph * (float)kBins);
        if (bin < 0) bin = 0;
        if (bin >= kBins) bin = kBins - 1;
        sum[bin] += d->cycleBuf[s];
        cnt[bin] += 1;
    }

    for (int i = 0; i < kBins; ++i) {
        d->newA[i] = (cnt[i] > 0) ? (sum[i] / (float)cnt[i]) : d->A[i];
    }

    if (!d->freeze) {
        float a = clampf(d->adaptAlpha, 0.0f, 1.0f);
        for (int i = 0; i < kBins; ++i) {
            d->A[i] += a * (d->newA[i] - d->A[i]);
        }
    }

    if (d->pendingCapture && d->locked) {
        for (int i = 0; i < kBins; ++i) d->A[i] = d->newA[i];
        d->freeze = true;
        d->pendingCapture = false;
    }

    combineS(d);
    d->herm.computeSlopes(d->S);

    for (int k = 0; k < kDisplayPoints; ++k) {
        float ph = (float)k / (float)(kDisplayPoints - 1);
        float pos = ph * (float)(cycleLen - 1);
        int i0 = (int)floorf(pos);
        int i1 = i0 + 1;
        if (i1 >= cycleLen) i1 = cycleLen - 1;
        float frac = pos - (float)i0;
        float xin = lerpf(d->cycleBuf[i0], d->cycleBuf[i1], frac);
        d->dispIn[k] = xin;

        float yLin = evalLinear(d->S, ph);
        float yCur = d->herm.eval(d->S, ph);
        d->dispOut[k] = lerpf(yLin, yCur, d->curve);
    }

    d->cycleWrite = 0;
}

// Simple LCG random, returns -1..1
static inline float lcgRand(uint32_t &state) {
    state = state * 1664525u + 1013904223u;
    return (float)(int32_t)state / (float)0x80000000u;
}

// Build sequencer target pattern into seqTarget[]
static void seqBuildPattern(struct _CessnaAwg_DTC *d) {
    switch (d->seqPattern) {
        case 0: // Ramp: 0→+1 across bins
            for (int i=0;i<kBins;i++)
                d->seqTarget[i] = ((float)i / (float)(kBins-1)) * 2.0f - 1.0f;
            break;
        case 1: // Triangle: 0→+1→0
            for (int i=0;i<kBins;i++) {
                float t = (float)i / (float)(kBins-1);
                d->seqTarget[i] = (t < 0.5f) ? (t*4.0f - 1.0f) : (3.0f - t*4.0f);
            }
            break;
        case 2: // Random: each bin gets an independent random walk step
            for (int i=0;i<kBins;i++) {
                d->seqTarget[i] = clampf(d->seqTarget[i] + lcgRand(d->seqRng) * d->seqDrift, -1.0f, 1.0f);
            }
            break;
        case 3: // Alternating: +1 -1 +1 -1...
            for (int i=0;i<kBins;i++)
                d->seqTarget[i] = (i & 1) ? -1.0f : 1.0f;
            break;
    }
    // Apply drift to non-random patterns too (random walk on top of the shape)
    if (d->seqPattern != 2 && d->seqDrift > 0.0f) {
        for (int i=0;i<kBins;i++) {
            d->seqTarget[i] = clampf(d->seqTarget[i] + lcgRand(d->seqRng) * d->seqDrift * 0.3f, -1.0f, 1.0f);
        }
    }
}

// Rotate the sequencer pattern by one bin (shift register advance)
static void seqRotate(struct _CessnaAwg_DTC *d) {
    float tmp = d->seqTarget[kBins-1];
    for (int i=kBins-1;i>0;i--) d->seqTarget[i] = d->seqTarget[i-1];
    d->seqTarget[0] = tmp;
}

// Called once per detected cycle to advance the sequencer
static void seqAdvance(struct _CessnaAwg_DTC *d) {
    d->seqCycleCount++;
    if (d->seqCycleCount < d->seqSpeed) return;
    d->seqCycleCount = 0;

    if (d->seqPattern == 2) {
        // Random: rebuild (random walk)
        seqBuildPattern(d);
    } else {
        // Others: rotate (shift register) + optional drift
        seqRotate(d);
        if (d->seqDrift > 0.0f) {
            for (int i=0;i<kBins;i++)
                d->seqTarget[i] = clampf(d->seqTarget[i] + lcgRand(d->seqRng) * d->seqDrift * 0.1f, -1.0f, 1.0f);
        }
    }
}

// Compute modX[] for one sample. For Ripple and Basis this is per-sample.
// For Sequencer, modX[] is only updated in seqAdvance (per cycle), but
// we still slew seqCurrent toward seqTarget here.
static __attribute__((noinline)) void computeModX(struct _CessnaAwg_DTC *d) {
    if (d->modMode == 0) {
        // Off — zero out (manual offsets still applied separately)
        for (int i=0;i<kBins;i++) d->modX[i] = 0.0f;
        return;
    }

    // Advance LFO
    d->lfoPhase += d->lfoPhaseInc;
    if (d->lfoPhase >= 1.0f) d->lfoPhase -= 1.0f;
    float lfoRad = d->lfoPhase * 2.0f * M_PI_F;

    if (d->modMode == 1) {
        // Ripple: each bin gets LFO with phase offset
        for (int i=0;i<kBins;i++) {
            float binPhase = lfoRad + d->modSpread * 2.0f * M_PI_F * ((float)i / (float)kBins);
            d->modX[i] = d->modDepth * sinf(binPhase);
        }
    } else if (d->modMode == 2) {
        // Basis: four shapes mixed
        // B1 = sine, B2 = cosine, B3 = saw/tilt, B4 = alternating
        // LFO drives coefficients with 90° phase offsets between them
        float lfoA = sinf(lfoRad);
        float lfoB = sinf(lfoRad + M_PI_F * 0.5f);
        float lfoC = sinf(lfoRad + M_PI_F);
        float lfoD = sinf(lfoRad + M_PI_F * 1.5f);

        for (int i=0;i<kBins;i++) {
            float t    = (float)i / (float)(kBins-1);
            float B1   = sinf(2.0f * M_PI_F * t);
            float B2   = cosf(2.0f * M_PI_F * t);
            float B3   = t * 2.0f - 1.0f;
            float B4   = (i & 1) ? -1.0f : 1.0f;

            // Static contribution
            float vstatic = d->basisStatic[0]*B1 + d->basisStatic[1]*B2
                          + d->basisStatic[2]*B3 + d->basisStatic[3]*B4;

            // LFO-driven contribution (coefficients modulated by LFO)
            float vlfo = lfoA*B1 + lfoB*B2 + lfoC*B3 + lfoD*B4;

            d->modX[i] = clampf(d->modDepth * lerpf(vstatic, vlfo, d->basisMix) * 0.25f, -1.0f, 1.0f);
        }
    } else {
        // Sequencer: slew seqCurrent toward seqTarget, output seqCurrent
        // Slew coefficient: smooth=0 → instant, smooth=1 → very slow
        float slew = 1.0f - d->seqSmooth * 0.999f;
        for (int i=0;i<kBins;i++) {
            d->seqCurrent[i] += slew * (d->seqTarget[i] - d->seqCurrent[i]);
            d->modX[i] = d->modDepth * d->seqCurrent[i];
        }
    }
}

// ----------------------------
// Sync mode helpers
// ----------------------------

// Shared correction: apply phase nudge and conditionally update frequency.
// Called by all three sync modes when they decide a correction should happen.
static inline void applySync(_CessnaAwg_DTC *d, float newPeriod) {
    // Phase correction (aggressiveness = lockStrength)
    if (!d->hardJump) {
        d->phase = d->phase * (1.0f - d->lockStrength);
    }
    // Frequency tracking (trackStrength controls how much phaseIncTarget updates)
    if (d->trackStrength > 0.0f) {
        float target = 1.0f / newPeriod;
        // Blend: full track → use new target directly; partial → interpolate
        d->phaseIncTarget = d->phaseIncTarget + d->trackStrength * (target - d->phaseIncTarget);
    }
    // (if trackStrength == 0, phaseIncTarget never changes → free-run)
}

// N-before-M: feed one lead/lag observation. Returns true if a correction fires.
// On firing, resets counters and calls applySync.
static bool nbmFeed(_CessnaAwg_DTC *d, bool isLead, float newPeriod) {
    if (isLead) d->nbmLeadCount++;
    else        d->nbmLagCount++;
    d->nbmTotalCount++;

    bool fired = false;

    // Condition 1: N of one type arrived before M total → output + reset
    if (d->nbmLeadCount >= d->nbmN || d->nbmLagCount >= d->nbmN) {
        applySync(d, newPeriod);
        fired = true;
        d->nbmLeadCount  = 0;
        d->nbmLagCount   = 0;
        d->nbmTotalCount = 0;
    }
    // Condition 2: M total arrived without N of either type → reset, no output
    else if (d->nbmTotalCount >= d->nbmM) {
        d->nbmLeadCount  = 0;
        d->nbmLagCount   = 0;
        d->nbmTotalCount = 0;
    }

    return fired;
}

// PLL: called once per detected edge.
// Separates frequency correction (PI loop on period error) from
// phase correction (gentle proportional nudge on phase error).
static void pllUpdate(_CessnaAwg_DTC *d, float newPeriod) {
    // --- Phase error ---
    // At the edge moment, d->phase should be ~0 if locked.
    // Wrap to -0.5..0.5: positive = oscillator ahead, negative = behind.
    float phErr = d->phase;
    if (phErr > 0.5f)  phErr -= 1.0f;
    if (phErr < -0.5f) phErr += 1.0f;

    // --- Frequency error ---
    // Difference between measured frequency and current estimate.
    float measuredFreq = 1.0f / newPeriod;
    float freqErr = measuredFreq - d->phaseIncTarget;

    // --- PI loop filter on frequency error ---
    // Kp scales pllBandwidth directly; Ki is much smaller to avoid windup.
    float Kp = d->pllBandwidth;
    float Ki = d->pllBandwidth * 0.05f;

    d->pllIntegrator += Ki * freqErr;
    // Clamp integrator to prevent windup (fixed: result now assigned back)
    if (d->pllIntegrator >  0.05f) d->pllIntegrator =  0.05f;
    if (d->pllIntegrator < -0.05f) d->pllIntegrator = -0.05f;

    float freqCorr = Kp * freqErr + d->pllIntegrator;

    // --- Update frequency estimate ---
    if (d->trackStrength > 0.0f) {
        float newTarget = d->phaseIncTarget + d->trackStrength * freqCorr;
        // Clamp to valid range
        float minInc = 1.0f / (float)d->maxPeriodSamps;
        float maxInc = 1.0f / (float)d->minPeriodSamps;
        if (newTarget < minInc) newTarget = minInc;
        if (newTarget > maxInc) newTarget = maxInc;
        d->phaseIncTarget = newTarget;
    }

    // --- Phase correction (separate from frequency) ---
    // Gentle proportional nudge — does not snap, just steers.
    // Scale by lockStrength and a conservative factor to avoid overcorrection.
    if (!d->hardJump) {
        float phCorr = d->lockStrength * phErr * 0.25f;
        d->phase -= phCorr;
        if (d->phase < 0.0f)  d->phase += 1.0f;
        if (d->phase >= 1.0f) d->phase -= 1.0f;
    }

    // Reset integrator on hard pitch jump to prevent windup on large transients
    if (d->hardJump) d->pllIntegrator = 0.0f;
}

static void step(_NT_algorithm *base, float *busFrames, int numFramesBy4) {
    auto *self = (_CessnaAwg*)base;
    if (!self || !self->dtc) return;
    auto *d = self->dtc;
    if (!self->v) return;

    const int numFrames = numFramesBy4 * 4;

    // Deferred grey-out — NT_algorithmIndex() is valid from step() onwards
    if (d->needsGrayOutUpdate) {
        d->needsGrayOutUpdate = false;
        updateModGrayOut(base);
    }

    const int inBus   = self->v[kParamInput];
    const int out1Bus = self->v[kParamOut1];
    const int out1Mode= self->v[kParamOut1Mode];
    const int out2Bus = self->v[kParamOut2];
    const int out2Mode= self->v[kParamOut2Mode];
    const int out3Bus = self->v[kParamOut3];
    const int out3Mode= self->v[kParamOut3Mode];

    float *in   = busPtr(busFrames, numFrames, inBus);
    float *out1 = busPtr(busFrames, numFrames, out1Bus);
    float *out2 = busPtr(busFrames, numFrames, out2Bus);
    float *out3 = busPtr(busFrames, numFrames, out3Bus);

    auto writeOut = [&](float* out, int n, float v, int mode) {
        if (!out) return;
        if (mode == 0) out[n] += v;
        else out[n] = v;
    };

    if (!in) {
        for (int n=0;n<numFrames;n++) {
            // Slew phaseInc toward target
            d->phaseInc += d->slewCoeff * (d->phaseIncTarget - d->phaseInc);

            computeModX(d);
            combineS(d);
            if (d->curve > 0.0001f) d->herm.computeSlopes(d->S);

            float lin = evalLinear(d->S, d->phase);
            float cur = (d->curve > 0.0001f) ? d->herm.eval(d->S, d->phase) : lin;
            float y = lerpf(lin, cur, d->curve);

            writeOut(out1, n, y, out1Mode);
            writeOut(out2, n, 0.0f, out2Mode);
            writeOut(out3, n, 0.0f, out3Mode);  // no input → error = 0

            d->phase += d->phaseInc;
            if (d->phase >= 1.0f) d->phase -= 1.0f;
        }
        return;
    }

    for (int n=0;n<numFrames;n++) {
        // Use raw bus value - NT internal float format, do not scale.
        float x = in[n];

        if (d->cycleWrite < kMaxCycleSamps) {
            d->cycleBuf[d->cycleWrite++] = x;
        }

        // Failsafe: if we never detect an edge and the buffer fills, reset state so we don't get stuck.
        if (d->cycleWrite >= kMaxCycleSamps) {
            d->cycleWrite = 0;
            d->wasBelowNeg = false;
            d->samplesSinceEdge = 0;
            d->locked = false;
            d->stableCount = 0;
        }

        bool edge = processCycleDetection(d, x);
        if (edge) {
            // cycleWrite was already incremented after writing x, so the current
            // edge sample is already in the buffer. Exclude it from the completed
            // cycle by using cycleWrite-1 as the length.
            int cycleLen = d->cycleWrite - 1;
            if (cycleLen < 0) cycleLen = 0;
            int p = d->lastPeriodSamples;
            if (p >= 8 && p <= kMaxCycleSamps) {
                if (cycleLen > p) cycleLen = p;
            }

            float newPeriod = d->smoothPeriod;

            if (d->syncMode == 0) {
                // Edge mode: correct immediately on every detected edge
                applySync(d, newPeriod);
            } else if (d->syncMode == 1) {
                // N-before-M: pre-condition with Edge behaviour until locked,
                // then switch to voting. This ensures the oscillator is roughly
                // aligned before N-before-M starts making correction decisions.
                if (!d->locked) {
                    // Not yet locked — use Edge to acquire
                    applySync(d, newPeriod);
                } else {
                    // Locked — feed N-before-M vote
                    bool isLead = (d->phase > 0.5f);
                    nbmFeed(d, isLead, newPeriod);
                }
            } else {
                // PLL: continuous frequency + gentle phase correction
                pllUpdate(d, newPeriod);
            }

            if (d->modMode == 3) seqAdvance(d);
            if (!d->hardJump) {
                finalizeCycle(d, cycleLen);
            }
            d->cycleWrite = 0;  // always reset buffer at edge
        }

        // Slew phaseInc toward target
        d->phaseInc += d->slewCoeff * (d->phaseIncTarget - d->phaseInc);

        computeModX(d);
        combineS(d);
        if (d->curve > 0.0001f) d->herm.computeSlopes(d->S);
        float yLin = evalLinear(d->S, d->phase);
        float yCur = (d->curve > 0.0001f) ? d->herm.eval(d->S, d->phase) : yLin;
        float y = lerpf(yLin, yCur, d->curve);

        writeOut(out1, n, y, out1Mode);
        writeOut(out2, n, x, out2Mode);
        writeOut(out3, n, x - y, out3Mode);  // error = input - resynth

        d->phase += d->phaseInc;
        if (d->phase >= 1.0f) d->phase -= 1.0f;
    }
}

// ----------------------------
// UI drawing
// ----------------------------
static bool draw(_NT_algorithm *base) {
    auto *self = (_CessnaAwg*)base;
    if (!self || !self->dtc) return false;
    auto *d = self->dtc;

    const int W = 256;
    const int H = 64;

    // --- Waveform area — below the standard parameter line (top 10px) ---
    int top    = 11;
    int bottom = H - 2;
    int plotH  = bottom - top;

    auto mapY = [&](float v)->int {
        float vv = clampf(v, -6.0f, 6.0f);
        float yn = 0.5f * (1.0f - (vv / 6.0f));
        int y = top + (int)roundf(yn * (float)plotH);
        if (y < top)    y = top;
        if (y >= H)     y = H-1;
        return y;
    };

    for (int n=1;n<kDisplayPoints;n++) {
        int x0 = (int)roundf(((float)(n-1) / (float)(kDisplayPoints-1)) * (float)(W-1));
        int x1 = (int)roundf(((float)n     / (float)(kDisplayPoints-1)) * (float)(W-1));

        int yIn0  = mapY(d->dispIn[n-1]);
        int yIn1  = mapY(d->dispIn[n]);
        int yOut0 = mapY(d->dispOut[n-1]);
        int yOut1 = mapY(d->dispOut[n]);

        NT_drawShapeI(kNT_line, x0, yIn0,  x1, yIn1,  6);
        NT_drawShapeI(kNT_line, x0, yOut0, x1, yOut1, 15);
    }

    // Lock indicator — small filled square, top-right of waveform area
    int lockCol = d->locked ? 15 : 3;
    NT_drawShapeI(kNT_rectangle, W-8, 11, W-2, 17, lockCol);

    return false;  // let NT draw the standard parameter line at top
}

static uint32_t hasCustomUi(_NT_algorithm *) { return 0; }
static void customUi(_NT_algorithm *, const _NT_uiData &) {}

static void calculateStaticRequirements(_NT_staticRequirements& req) { req.dram = 0; }

static void initialise(_NT_staticMemoryPtrs& ptrs, const _NT_staticRequirements& req) {
    (void)ptrs; (void)req;
}

// IMPORTANT: calculateRequirements is already defined earlier in the file.
// Do NOT redefine it here.

static _NT_factory g_factory = {
    .guid = NT_MULTICHAR('C','s','n','A'),
    .name = "Cessna++ AWG",
    .description = "Phase-segmented, servo-tracking arbitrary waveform generator",
    .numSpecifications = 0,
    .specifications = nullptr,

    .calculateStaticRequirements = calculateStaticRequirements,
    .initialise = initialise,
    .calculateRequirements = calculateRequirements,
    .construct = construct,
    .parameterChanged = parameterChanged,
    .step = step,
    .draw = draw,
    .midiRealtime = nullptr,
    .midiMessage = nullptr,
    .tags = kNT_tagUtility,
    .hasCustomUi = hasCustomUi,
    .customUi = customUi,
    .setupUi = nullptr,
    .serialise = nullptr,
    .deserialise = nullptr,
    .midiSysEx = nullptr,
    .parameterUiPrefix = nullptr,
    .parameterString = nullptr,
};

extern "C" uintptr_t pluginEntry(_NT_selector selector, uint32_t data) {
    (void)data;
    switch (selector) {
        case kNT_selector_version:      return (uintptr_t)kNT_apiVersionCurrent;
        case kNT_selector_numFactories: return 1;
        case kNT_selector_factoryInfo:  return (uintptr_t)&g_factory;
        default:                        return 0;
    }
}