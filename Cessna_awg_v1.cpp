/*
 * Cessna++ AWG (v1) – Phase-Segmented, Servo-Tracking Arbitrary Waveform Generator
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
 *
 * Notes:
 * - This is a first-pass implementation skeleton intended to be refined against the official API
 *   examples and build system.
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

// Auto-threshold parameters
static constexpr float kEnvBeta        = 0.0015f;  // envelope smoothing (per sample)
static constexpr float kThreshK        = 0.10f;    // threshold = max(Tmin, k*env)
static constexpr float kThreshMin      = 0.0025f;  // minimum threshold in input sample units

// Tick highlight
static constexpr float kTickChangeThresh = 0.01f;   // change threshold in X[i] (assuming [-1,1])
static constexpr float kTickDecayPerDraw = 0.82f;   // per draw call decay factor

// Max samples to keep for one captured input cycle
// (at 48kHz and 20Hz minimum, one cycle is 2400 samples; add a little margin)
static constexpr int kMaxCycleSamps = (int)(kDefaultSampleRate / kMinHz) + 64;

// ----------------------------
// Helpers
// ----------------------------
static inline float clampf(float x, float lo, float hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

static inline float lerpf(float a, float b, float t) {
    return a + (b - a) * t;
}

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
// This uses a Fritsch–Carlson style slope limiting for monotone segments.
// For segments that are not monotone, it tends toward Catmull slopes but clamps overshoot.
struct PeriodicHermite16 {
    float m[kBins] = {0}; // slopes at points

    void computeSlopes(const float *p) {
        // Finite differences (periodic)
        float d[kBins];
        for (int i=0;i<kBins;i++) {
            float p0 = p[i];
            float p1 = p[(i+1)&15];
            d[i] = (p1 - p0);
        }

        // Initial slope guess: average adjacent deltas (Catmull-ish)
        for (int i=0;i<kBins;i++) {
            float dm = d[(i-1)&15];
            float dp = d[i];
            m[i] = 0.5f * (dm + dp);
        }

        // Monotone limiting per segment
        for (int i=0;i<kBins;i++) {
            float dp = d[i];
            float a = m[i];
            float b = m[(i+1)&15];

            if (dp == 0.0f) {
                // Flat segment -> slopes should be zero to avoid wiggle
                m[i] = 0.0f;
                m[(i+1)&15] = 0.0f;
                continue;
            }

            // If segment is monotone, enforce Fritsch–Carlson style constraint
            // to prevent overshoot: a/d and b/d must be nonnegative and (a/d)^2+(b/d)^2 <= 9
            float r1 = a / dp;
            float r2 = b / dp;

            if (r1 < 0.0f) m[i] = 0.0f;
            if (r2 < 0.0f) m[(i+1)&15] = 0.0f;

            r1 = m[i] / dp;
            r2 = m[(i+1)&15] / dp;

            float sumsq = r1*r1 + r2*r2;
            if (sumsq > 9.0f) {
                float t = 3.0f / sqrtf(sumsq);
                m[i] *= t;
                m[(i+1)&15] *= t;
            }
        }
    }

    inline float eval(const float *p, float phase01) const {
        // phase01 in [0,1)
        float x = phase01 * (float)kBins;
        int i = (int)floorf(x);
        float t = x - (float)i;
        i &= 15;
        int i1 = (i+1) & 15;

        float p0 = p[i];
        float p1 = p[i1];
        float m0 = m[i];
        float m1 = m[i1];

        // Hermite basis
        float t2 = t*t;
        float t3 = t2*t;
        float h00 =  2.0f*t3 - 3.0f*t2 + 1.0f;
        float h10 =        t3 - 2.0f*t2 + t;
        float h01 = -2.0f*t3 + 3.0f*t2;
        float h11 =        t3 -       t2;

        // For unit spacing between control points
        return h00*p0 + h10*m0 + h01*p1 + h11*m1;
    }
};

// ----------------------------
// DTC (fast memory)
// ----------------------------
struct _CessnaAwg_DTC {
    // Control / state
    float sampleRate;

    // Parameters
    float depthD;          // global deviation depth
    float adaptAlpha;      // servo convergence (0..1)
    float curve;           // 0=linear, 1=curved
    float lockStrength;    // 0..1 for future phase-locking nuance; 1 implies hard reset
    bool  freeze;
    bool  pendingCapture;  // user requested capture (armed)

    // I2C-driven deviations (mapped to parameters)
    float X[kBins];
    float X_prev[kBins];

    // Tick highlight
    float tickHi[kBins];

    // Cycle detection
    DcBlocker dc;
    float env;
    float thresh;
    bool  wasBelowNeg;
    int   samplesSinceEdge;
    int   lastPeriodSamples;
    float smoothPeriod;
    bool  locked;

    int   prevPeriodSamples;
    int   stableCount;

    // One-cycle capture buffer for analysis + display
    float cycleBuf[kMaxCycleSamps];
    int   cycleWrite;
    int   minPeriodSamps;
    int   maxPeriodSamps;

    // Stored analyzed waveform bins
    float A[kBins];

    // Temporary new bins (computed per cycle)
    float newA[kBins];

    // Resynth bins after deviation
    float S[kBins];

    // Hermite slopes for curved interpolation
    PeriodicHermite16 herm;

    // Oscillator phase
    float phase;           // [0,1)
    float phaseInc;        // per-sample

    // Display buffers (phase-normalised)
    float dispIn[kDisplayPoints];
    float dispOut[kDisplayPoints];
    bool  dispInit;

    _CessnaAwg_DTC() {
        memset(this, 0, sizeof(*this));
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

    kParamDepth,
    kParamAdapt,
    kParamCurve,
    kParamLock,
    kParamFreeze,
    kParamCapture,

    // 16 deviation params (to be mapped from I2C)
    kParamX1,
    // ... sequential
    kParamX16 = kParamX1 + 15,
};

static const char* onOffStrings[] = {"Off", "On", nullptr};

// ----------------------------
// Requirements / construct
// ----------------------------
static constexpr int kNumParams = (kParamX16 + 1);

// ----------------------------
// Parameters (new API)
// ----------------------------
static const _NT_parameter g_parameters[kNumParams] = {
    // IO
    NT_PARAMETER_AUDIO_INPUT("in", 0, 1)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("out1", 0, 13)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("out2", 0, 14)
    // Main controls
    { .name="Depth",   .min=0,   .max=100, .def=0,   .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Adapt",   .min=0,   .max=100, .def=100, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Curve",   .min=0,   .max=100, .def=0,   .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Lock",    .min=0,   .max=100, .def=100, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="freeze",  .min=0,   .max=1,   .def=0,   .unit=kNT_unitEnum,    .scaling=kNT_scalingNone, .enumStrings=onOffStrings },
    { .name="capture", .min=0,   .max=1,   .def=0,   .unit=kNT_unitEnum,    .scaling=kNT_scalingNone, .enumStrings=onOffStrings },

    // 16 deviation controls (percent)
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
};

static const uint8_t g_pageRoutingIdx[] = {
    (uint8_t)kParamInput,
    (uint8_t)kParamOut1, (uint8_t)kParamOut1Mode,
    (uint8_t)kParamOut2, (uint8_t)kParamOut2Mode,
};

static const uint8_t g_pageMainIdx[] = {
    (uint8_t)kParamDepth,
    (uint8_t)kParamAdapt,
    (uint8_t)kParamCurve,
    (uint8_t)kParamLock,
    (uint8_t)kParamFreeze,
    (uint8_t)kParamCapture,
};

static const uint8_t g_pageXIdx[] = {
    (uint8_t)kParamX1,  (uint8_t)(kParamX1+1),  (uint8_t)(kParamX1+2),  (uint8_t)(kParamX1+3),
    (uint8_t)(kParamX1+4),  (uint8_t)(kParamX1+5),  (uint8_t)(kParamX1+6),  (uint8_t)(kParamX1+7),
    (uint8_t)(kParamX1+8),  (uint8_t)(kParamX1+9),  (uint8_t)(kParamX1+10), (uint8_t)(kParamX1+11),
    (uint8_t)(kParamX1+12), (uint8_t)(kParamX1+13), (uint8_t)(kParamX1+14), (uint8_t)(kParamX1+15),
};

static const _NT_parameterPage g_pages_arr[] = {
    { .name="Routing", .numParams=(uint8_t)ARRAY_SIZE(g_pageRoutingIdx), .group=0, .unused={0,0}, .params=g_pageRoutingIdx },
    { .name="Main",    .numParams=(uint8_t)ARRAY_SIZE(g_pageMainIdx),    .group=0, .unused={0,0}, .params=g_pageMainIdx },
    { .name="Offsets", .numParams=(uint8_t)ARRAY_SIZE(g_pageXIdx),       .group=1, .unused={0,0}, .params=g_pageXIdx },
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



    d->sampleRate = (NT_globals.sampleRate > 0) ? (float)NT_globals.sampleRate : kDefaultSampleRate;

    d->depthD = 0.0f;
    d->adaptAlpha = 1.0f;
    d->curve = 0.0f;
    d->lockStrength = 1.0f;
    d->freeze = false;
    d->pendingCapture = false;

    d->env = 0.0f;
    d->thresh = kThreshMin;
    d->wasBelowNeg = false;
    d->samplesSinceEdge = 0;
    d->lastPeriodSamples = (int)(d->sampleRate / 110.0f);
    d->smoothPeriod = (float)d->lastPeriodSamples;
    d->locked = false;
    d->prevPeriodSamples = d->lastPeriodSamples;
    d->stableCount = 0;

    d->minPeriodSamps = (int)floorf(d->sampleRate / kMaxHz);
    d->maxPeriodSamps = (int)ceilf(d->sampleRate / kMinHz);
    if (d->minPeriodSamps < 4) d->minPeriodSamps = 4;

    // Init A as a sine-ish default
    for (int i=0;i<kBins;i++) {
        float ph = ((float)i + 0.5f) / (float)kBins;
        d->A[i] = sinf(2.0f*M_PI_F*ph);
        d->X[i] = 0.0f;
        d->X_prev[i] = 0.0f;
        d->tickHi[i] = 0.0f;
    }

    // Clear cycle buffer
    for (int i=0;i<kMaxCycleSamps;i++){ d->cycleBuf[i]=0.0f; }
    d->cycleWrite = 0;

    d->phase = 0.0f;
    d->phaseInc = 1.0f / d->smoothPeriod;

    // Display init
    for (int i=0;i<kDisplayPoints;i++){ d->dispIn[i]=0.0f; d->dispOut[i]=0.0f; }
    // Seed the resynth display with the initial sine-ish bins.
    for (int i=0;i<kBins;i++) d->S[i] = d->A[i];
    d->herm.computeSlopes(d->S);
    for (int k=0;k<kDisplayPoints;k++) {
        float ph = (float)k / (float)(kDisplayPoints - 1);
        float yLin = evalLinear(d->S, ph);
        float yCur = d->herm.eval(d->S, ph);
        d->dispOut[k] = lerpf(yLin, yCur, d->curve);
    }
    d->dispInit = true;

    return self;
}

// ----------------------------
// Parameter changes
// ----------------------------

// ----------------------------
// Parameter changes (new API)
// ----------------------------
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

    switch (p) {
        case kParamDepth:   d->depthD      = getPct01(kParamDepth); break;
        case kParamAdapt:   d->adaptAlpha  = getPct01(kParamAdapt); break;
        case kParamCurve:   d->curve       = getPct01(kParamCurve); break;
        case kParamLock:    d->lockStrength= getPct01(kParamLock); break;
        case kParamFreeze:  d->freeze      = (self->v[kParamFreeze] != 0); break;
        case kParamCapture: d->pendingCapture = (self->v[kParamCapture] != 0); break;
        default:
            if (p >= kParamX1 && p <= kParamX16) {
                int i = p - kParamX1;
                d->X[i] = getPctSigned(p);
            }
            break;
    }
}


static inline float* busPtr(float* busFrames, int numFrames, int busParamValue) {
    // busParamValue: 0 = none, otherwise 1..kNT_lastBus
    if (!busFrames) return nullptr;
    if (busParamValue <= 0) return nullptr;
    int busIndex = busParamValue - 1;
    if (busIndex < 0 || busIndex >= 28) return nullptr;
    return busFrames + (busIndex * numFrames);
}

// ----------------------------
// Audio step (new API)
// ----------------------------
// Forward declarations (helpers defined later)
static inline float evalLinear(const float *bins, float phase);
static bool processCycleDetection(struct _CessnaAwg_DTC *d, float x);
static void finalizeCycle(struct _CessnaAwg_DTC *d, int cycleLen);

// ----------------------------
// Helper implementations
// ----------------------------

// Linear interpolation over 16 periodic bins. phase in [0,1).
static inline float evalLinear(const float *bins, float phase) {
    // Keep phase in [0,1)
    phase -= floorf(phase);
    float x = phase * (float)kBins;
    int i = (int)floorf(x);
    float f = x - (float)i;
    int i0 = i & (kBins - 1);
    int i1 = (i + 1) & (kBins - 1);
    return lerpf(bins[i0], bins[i1], f);
}

// Rising zero-crossing with auto hysteresis threshold.
// Returns true once per detected cycle boundary.
static bool processCycleDetection(struct _CessnaAwg_DTC *d, float x) {
    // DC-block for more reliable crossings
    float y = d->dc.process(x);

    // Envelope follower (abs)
    float ay = fabsf(y);
    d->env = (1.0f - kEnvBeta) * d->env + kEnvBeta * ay;

    // Auto threshold
    d->thresh = fmaxf(kThreshMin, kThreshK * d->env);

    // Hysteresis: require we went sufficiently negative before accepting a positive crossing
    if (y < -d->thresh) d->wasBelowNeg = true;

    d->samplesSinceEdge++;

    if (d->wasBelowNeg && (y > d->thresh)) {
        // Candidate rising edge
        int period = d->samplesSinceEdge;
        d->samplesSinceEdge = 0;
        d->wasBelowNeg = false;

        // Validate plausible frequency
        if (period >= d->minPeriodSamps && period <= d->maxPeriodSamps) {
            d->lastPeriodSamples = period;

            // Smooth the estimate a bit to avoid jitter
            const float beta = 0.10f;
            d->smoothPeriod = (1.0f - beta) * d->smoothPeriod + beta * (float)period;
            if (d->smoothPeriod < 4.0f) d->smoothPeriod = 4.0f;
            d->phaseInc = 1.0f / d->smoothPeriod;

            // Lock heuristic: period stability over several consecutive cycles
            float rel = fabsf((float)period - (float)d->prevPeriodSamples) / (float)d->prevPeriodSamples;
            if (rel < 0.02f) d->stableCount++; else d->stableCount = 0;
            d->prevPeriodSamples = period;
            d->locked = (d->stableCount >= 2);

            return true;
        }
    }
    return false;
}

// Finalize a captured cycle in d->cycleBuf[0..cycleLen): compute per-bin means,
// update A with servo (unless frozen), optional capture->freeze, and refresh display.
static void finalizeCycle(struct _CessnaAwg_DTC *d, int cycleLen) {
    if (cycleLen < 8) {
        d->cycleWrite = 0;
        return;
    }

    // Compute new bin means from the actual captured cycle length.
    float sum[kBins] = {0};
    int   cnt[kBins] = {0};

    for (int s = 0; s < cycleLen; ++s) {
        float ph = (float)s / (float)cycleLen; // [0,1)
        int bin = (int)floorf(ph * (float)kBins);
        if (bin < 0) bin = 0;
        if (bin >= kBins) bin = kBins - 1;
        sum[bin] += d->cycleBuf[s];
        cnt[bin] += 1;
    }

    for (int i = 0; i < kBins; ++i) {
        d->newA[i] = (cnt[i] > 0) ? (sum[i] / (float)cnt[i]) : d->A[i];
    }

    // Servo update of stored bins unless frozen
    if (!d->freeze) {
        float a = clampf(d->adaptAlpha, 0.0f, 1.0f);
        for (int i = 0; i < kBins; ++i) {
            d->A[i] += a * (d->newA[i] - d->A[i]);
        }
    }

    // Capture only when locked: overwrite A with the freshly analysed cycle, then Freeze ON
    if (d->pendingCapture && d->locked) {
        for (int i = 0; i < kBins; ++i) d->A[i] = d->newA[i];
        d->freeze = true;
        d->pendingCapture = false;
    }

    // Update S and slopes for display render
    float D = d->depthD;
    for (int i = 0; i < kBins; ++i) d->S[i] = d->A[i] + D * d->X[i];
    d->herm.computeSlopes(d->S);

    // Refresh display buffers: resample the captured input cycle to fixed width.
    for (int k = 0; k < kDisplayPoints; ++k) {
        float ph = (float)k / (float)(kDisplayPoints - 1);
        float pos = ph * (float)(cycleLen - 1);
        int i0 = (int)floorf(pos);
        int i1 = i0 + 1;
        if (i1 >= cycleLen) i1 = cycleLen - 1;
        float frac = pos - (float)i0;
        float xin = lerpf(d->cycleBuf[i0], d->cycleBuf[i1], frac);
        d->dispIn[k] = xin;

        // Resynth one cycle directly from the current S[] shape.
        float yLin = evalLinear(d->S, ph);
        float yCur = d->herm.eval(d->S, ph);
        d->dispOut[k] = lerpf(yLin, yCur, d->curve);
    }

    d->cycleWrite = 0;
}

static void step(_NT_algorithm *base, float *busFrames, int numFramesBy4) {
    auto *self = (_CessnaAwg*)base;
    if (!self || !self->dtc) return;
    auto *d = self->dtc;
    if (!self->v) return;

    const int numFrames = numFramesBy4 * 4;

    const int inBus  = self->v[kParamInput];
    const int out1Bus = self->v[kParamOut1];
    const int out1Mode = self->v[kParamOut1Mode];
    const int out2Bus = self->v[kParamOut2];
    const int out2Mode = self->v[kParamOut2Mode];

    float *in   = busPtr(busFrames, numFrames, inBus);
    float *out1 = busPtr(busFrames, numFrames, out1Bus);
    float *out2 = busPtr(busFrames, numFrames, out2Bus);

    auto writeOut = [&](float* out, int n, float v, int mode) {
        if (!out) return;
        // Convention: 0 = Add, 1 = Replace (matches many built-in algorithms)
        if (mode == 0) out[n] += v;
        else out[n] = v;
    };

    if (!in) {
        // No input: just output resynth based on current bins
        for (int n=0;n<numFrames;n++) {
            float D = d->depthD;
            for (int i=0;i<kBins;i++) d->S[i] = d->A[i] + D*d->X[i];
            d->herm.computeSlopes(d->S);

            float lin = evalLinear(d->S, d->phase);
            float cur = d->herm.eval(d->S, d->phase);
            float y = lerpf(lin, cur, d->curve);

            writeOut(out1, n, y, out1Mode);
            writeOut(out2, n, 0.0f, out2Mode);

            d->phase += d->phaseInc;
            if (d->phase >= 1.0f) d->phase -= 1.0f;
        }
        return;
    }

    for (int n=0;n<numFrames;n++) {
        float x = in[n];

        // Tick highlight on change of offsets (cheap and fine at audio-rate for 16 params)
        for (int i=0;i<kBins;i++) {
            float dx = fabsf(d->X[i] - d->X_prev[i]);
            if (dx > kTickChangeThresh) d->tickHi[i] = 1.0f;
            d->X_prev[i] = d->X[i];
        }

        // Capture input samples into the cycle buffer (bounded)
        if (d->cycleWrite < kMaxCycleSamps) {
            d->cycleBuf[d->cycleWrite++] = x;
        }

        // Cycle detection
        bool edge = processCycleDetection(d, x);
        if (edge) {
            int cycleLen = d->cycleWrite;
            // Prefer the measured period if it matches the buffer length reasonably.
            int p = d->lastPeriodSamples;
            if (p >= 8 && p <= kMaxCycleSamps) {
                // If the buffer overran slightly (or underfilled), clamp to measured period.
                if (cycleLen > p) cycleLen = p;
                // If severely underfilled, just use what we have.
            }

            // If lock strength == 100%, hard reset phase at cycle boundary
            if (d->lockStrength >= 0.999f) d->phase = 0.0f;

            finalizeCycle(d, cycleLen);
        }

        // Build S[] for audio-rate resynth (so offsets respond immediately).
        float D = d->depthD;
        for (int i=0;i<kBins;i++) d->S[i] = d->A[i] + D * d->X[i];
        // Only compute slopes if we need curved interpolation.
        if (d->curve > 0.0001f) d->herm.computeSlopes(d->S);

        float yLin = evalLinear(d->S, d->phase);
        float yCur = (d->curve > 0.0001f) ? d->herm.eval(d->S, d->phase) : yLin;
        float y = lerpf(yLin, yCur, d->curve);

        writeOut(out1, n, y, out1Mode);
        writeOut(out2, n, x, out2Mode); // Dry

        // Advance phase
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

    const int width = 256;
    const int height = 64;

    // Decay tick highlights
    for (int i=0;i<kBins;i++) d->tickHi[i] *= kTickDecayPerDraw;

    // Draw baseline
    NT_drawShapeI(kNT_line, 0, height-1, width-1, height-1, 2);

    // Draw bottom ticks + labels 1..16 at bin centers
    for (int i=0;i<kBins;i++) {
        float xc = ((float)i + 0.5f) / (float)kBins;
        int x = (int)roundf(xc * (float)(width-1));
        int y1 = height - 1;
        int y0 = height - 8;
        int col = (int)lerpf(4.0f, 15.0f, clampf(d->tickHi[i], 0.0f, 1.0f));
        NT_drawShapeI(kNT_line, x, y0, x, y1, col);

        // Tiny label
        char buf[4];
        int label = i + 1;
        if (label < 10) {
            buf[0] = '0' + label;
            buf[1] = '\0';
        } else {
            buf[0] = '1';
            buf[1] = '0' + (label - 10);
            buf[2] = '\0';
        }
        // Draw text slightly above bottom; API provides NT_drawText? Use NT_drawTextI if present.
        // Fallback: omit if unavailable.
#ifdef NT_drawTextI
        NT_drawTextI(x-2, height-16, buf, col);
#endif
    }

    // Draw input and resynth waveforms (overlay)
    // Map amplitude [-1,1] to y [0, height-20]
    int top = 2;
    int bottom = height - 18;
    int plotH = bottom - top;

    auto mapY = [&](float v)->int {
        float vv = clampf(v, -1.2f, 1.2f);
        float yn = 0.5f * (1.0f - (vv / 1.2f));
        int y = top + (int)roundf(yn * (float)plotH);
        if (y < 0) y = 0;
        if (y >= height) y = height-1;
        return y;
    };

    // Draw polylines
    for (int n=1;n<kDisplayPoints;n++) {
        int x0 = (int)roundf(((float)(n-1) / (float)(kDisplayPoints-1)) * (float)(width-1));
        int x1 = (int)roundf(((float)n     / (float)(kDisplayPoints-1)) * (float)(width-1));

        int yIn0  = mapY(d->dispIn[n-1]);
        int yIn1  = mapY(d->dispIn[n]);
        int yOut0 = mapY(d->dispOut[n-1]);
        int yOut1 = mapY(d->dispOut[n]);

        NT_drawShapeI(kNT_line, x0, yIn0,  x1, yIn1,  6);  // input (dim)
        NT_drawShapeI(kNT_line, x0, yOut0, x1, yOut1, 15); // resynth (bright)
    }

    // Optional: lock indicator (small dot)
    int lockCol = d->locked ? 15 : 4;
    NT_drawShapeI(kNT_rectangle, width-8, 2, width-2, 8, lockCol);

    return true;
}

static uint32_t hasCustomUi(_NT_algorithm *) {
    return 0; // no custom controls (yet)
}

static void customUi(_NT_algorithm *, const _NT_uiData &) {}


// ----------------------------
// Factory + plug-in entry (API v4+ style)
// ----------------------------

static void calculateStaticRequirements(_NT_staticRequirements& req) {
    // No shared allocations.
    req.dram = 0;
}

static void initialise(_NT_staticMemoryPtrs& ptrs, const _NT_staticRequirements& req) {
    (void)ptrs;
    (void)req;
}

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
