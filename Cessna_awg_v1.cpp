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
static constexpr float kEnvBeta        = 0.0015f;  // envelope smoothing (per sample)

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

// ---- Forward declarations (MUST appear before construct/step) ----
static inline float evalLinear(const float *bins, float phase01);
static bool processCycleDetection(struct _CessnaAwg_DTC *d, float x);
static void finalizeCycle(struct _CessnaAwg_DTC *d, int cycleLen);

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

        for (int i=0;i<kBins;i++) {
            float dp = d[i];
            float a = m[i];
            float b = m[(i+1)&15];

            if (dp == 0.0f) {
                m[i] = 0.0f;
                m[(i+1)&15] = 0.0f;
                continue;
            }

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
    bool  freeze;
    bool  pendingCapture;

    float X[kBins];
    float X_prev[kBins];

    float tickHi[kBins];

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

    float dispIn[kDisplayPoints];
    float dispOut[kDisplayPoints];
    bool  dispInit;

    _CessnaAwg_DTC() {
        sampleRate = kDefaultSampleRate;

        depthD = 0.0f;
        adaptAlpha = 1.0f;
        curve = 0.0f;
        lockStrength = 1.0f;
        freeze = false;
        pendingCapture = false;

        env = 0.0f;
        thresh = 0.0f;
        wasBelowNeg = false;
        samplesSinceEdge = 0;
        lastPeriodSamples = (int)(kDefaultSampleRate / 110.0f);
        smoothPeriod = (float)lastPeriodSamples;
        locked = false;

        prevPeriodSamples = lastPeriodSamples;
        stableCount = 0;

        cycleWrite = 0;
        minPeriodSamps = (int)floorf(sampleRate / kMaxHz);
        maxPeriodSamps = (int)ceilf(sampleRate / kMinHz);
        if (minPeriodSamps < 4) minPeriodSamps = 4;

        phase = 0.0f;
        phaseInc = 1.0f / smoothPeriod;

        dispInit = false;

        for (int i=0;i<kBins;i++) {
            X[i] = 0.0f;
            X_prev[i] = 0.0f;
            tickHi[i] = 0.0f;

            float ph = ((float)i + 0.5f) / (float)kBins;
            A[i] = sinf(2.0f*M_PI_F*ph);
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

    kParamDepth,
    kParamAdapt,
    kParamCurve,
    kParamLock,
    kParamFreeze,
    kParamCapture,

    kParamX1,
    kParamX16 = kParamX1 + 15,
};

static const char* onOffStrings[] = {"Off", "On", nullptr};
static constexpr int kNumParams = (kParamX16 + 1);

// ----------------------------
// Parameters (new API)
// ----------------------------
static const _NT_parameter g_parameters[kNumParams] = {
    NT_PARAMETER_AUDIO_INPUT("in", 0, 1)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("out1", 0, 13)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("out2", 0, 14)

    { .name="Depth",   .min=0,   .max=100, .def=0,   .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Adapt",   .min=0,   .max=100, .def=100, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Curve",   .min=0,   .max=100, .def=0,   .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Lock",    .min=0,   .max=100, .def=100, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="freeze",  .min=0,   .max=1,   .def=0,   .unit=kNT_unitEnum,    .scaling=kNT_scalingNone, .enumStrings=onOffStrings },
    { .name="capture", .min=0,   .max=1,   .def=0,   .unit=kNT_unitEnum,    .scaling=kNT_scalingNone, .enumStrings=onOffStrings },

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
    d->minPeriodSamps = (int)floorf(d->sampleRate / kMaxHz);
    d->maxPeriodSamps = (int)ceilf(d->sampleRate / kMinHz);
    if (d->minPeriodSamps < 4) d->minPeriodSamps = 4;

    d->smoothPeriod = (float)d->lastPeriodSamples;
    d->phaseInc = 1.0f / d->smoothPeriod;

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
        case kParamDepth:   d->depthD       = getPct01(kParamDepth); break;
        case kParamAdapt:   d->adaptAlpha   = getPct01(kParamAdapt); break;
        case kParamCurve:   d->curve        = getPct01(kParamCurve); break;
        case kParamLock:    d->lockStrength = getPct01(kParamLock); break;
        case kParamFreeze:  d->freeze       = (self->v[kParamFreeze] != 0); break;
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
 * Revised cycle detector:
 * - Use DC-blocked signal ONLY for envelope estimation.
 * - Use RAW input x for hysteresis crossing (much more robust).
 * - More permissive threshold: thresh = max(0.0005, 0.05 * env)
 * - Loosen lock criterion (TEMP): 5% period tolerance, 1 stable cycle.
 */
static bool processCycleDetection(struct _CessnaAwg_DTC *d, float x) {
    // DC-block only for envelope estimate (stable against DC offsets)
    float yEnv = d->dc.process(x);

    float ay = fabsf(yEnv);
    d->env = (1.0f - kEnvBeta) * d->env + kEnvBeta * ay;

    // More permissive auto-threshold (still adapts to level)
    const float thresh = fmaxf(0.0005f, 0.05f * d->env);
    d->thresh = thresh;

    // Use RAW x for crossing/hysteresis
    if (x < -thresh) d->wasBelowNeg = true;

    d->samplesSinceEdge++;

    if (d->wasBelowNeg && (x > thresh)) {
        int period = d->samplesSinceEdge;
        d->samplesSinceEdge = 0;
        d->wasBelowNeg = false;

        if (period >= d->minPeriodSamps && period <= d->maxPeriodSamps) {
            d->lastPeriodSamples = period;

            const float beta = 0.10f;
            d->smoothPeriod = (1.0f - beta) * d->smoothPeriod + beta * (float)period;
            if (d->smoothPeriod < 4.0f) d->smoothPeriod = 4.0f;
            d->phaseInc = 1.0f / d->smoothPeriod;

            // TEMP lock heuristic: easier to reach while we validate on-device
            float rel = fabsf((float)period - (float)d->prevPeriodSamples) / (float)d->prevPeriodSamples;
            if (rel < 0.05f) d->stableCount++; else d->stableCount = 0;
            d->prevPeriodSamples = period;
            d->locked = (d->stableCount >= 1);

            return true;
        }
    }
    return false;
}

static void finalizeCycle(struct _CessnaAwg_DTC *d, int cycleLen) {
    if (cycleLen < 8) {
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

    float D = d->depthD;
    for (int i = 0; i < kBins; ++i) d->S[i] = d->A[i] + D * d->X[i];
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

static void step(_NT_algorithm *base, float *busFrames, int numFramesBy4) {
    auto *self = (_CessnaAwg*)base;
    if (!self || !self->dtc) return;
    auto *d = self->dtc;
    if (!self->v) return;

    const int numFrames = numFramesBy4 * 4;

    const int inBus   = self->v[kParamInput];
    const int out1Bus = self->v[kParamOut1];
    const int out1Mode= self->v[kParamOut1Mode];
    const int out2Bus = self->v[kParamOut2];
    const int out2Mode= self->v[kParamOut2Mode];

    float *in   = busPtr(busFrames, numFrames, inBus);
    float *out1 = busPtr(busFrames, numFrames, out1Bus);
    float *out2 = busPtr(busFrames, numFrames, out2Bus);

    auto writeOut = [&](float* out, int n, float v, int mode) {
        if (!out) return;
        if (mode == 0) out[n] += v;
        else out[n] = v;
    };

    if (!in) {
        for (int n=0;n<numFrames;n++) {
            float D = d->depthD;
            for (int i=0;i<kBins;i++) d->S[i] = d->A[i] + D*d->X[i];
            if (d->curve > 0.0001f) d->herm.computeSlopes(d->S);

            float lin = evalLinear(d->S, d->phase);
            float cur = (d->curve > 0.0001f) ? d->herm.eval(d->S, d->phase) : lin;
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

        for (int i=0;i<kBins;i++) {
            float dx = fabsf(d->X[i] - d->X_prev[i]);
            if (dx > kTickChangeThresh) d->tickHi[i] = 1.0f;
            d->X_prev[i] = d->X[i];
        }

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
            int cycleLen = d->cycleWrite;
            int p = d->lastPeriodSamples;
            if (p >= 8 && p <= kMaxCycleSamps) {
                if (cycleLen > p) cycleLen = p;
            }

            if (d->lockStrength >= 0.999f) d->phase = 0.0f;
            finalizeCycle(d, cycleLen);
        }

        float D = d->depthD;
        for (int i=0;i<kBins;i++) d->S[i] = d->A[i] + D * d->X[i];
        if (d->curve > 0.0001f) d->herm.computeSlopes(d->S);

        float yLin = evalLinear(d->S, d->phase);
        float yCur = (d->curve > 0.0001f) ? d->herm.eval(d->S, d->phase) : yLin;
        float y = lerpf(yLin, yCur, d->curve);

        writeOut(out1, n, y, out1Mode);
        writeOut(out2, n, x, out2Mode);

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

    for (int i=0;i<kBins;i++) d->tickHi[i] *= kTickDecayPerDraw;

    NT_drawShapeI(kNT_line, 0, height-1, width-1, height-1, 2);

    for (int i=0;i<kBins;i++) {
        float xc = ((float)i + 0.5f) / (float)kBins;
        int x = (int)roundf(xc * (float)(width-1));
        int y1 = height - 1;
        int y0 = height - 8;
        int col = (int)lerpf(4.0f, 15.0f, clampf(d->tickHi[i], 0.0f, 1.0f));
        NT_drawShapeI(kNT_line, x, y0, x, y1, col);
    }

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

    for (int n=1;n<kDisplayPoints;n++) {
        int x0 = (int)roundf(((float)(n-1) / (float)(kDisplayPoints-1)) * (float)(width-1));
        int x1 = (int)roundf(((float)n     / (float)(kDisplayPoints-1)) * (float)(width-1));

        int yIn0  = mapY(d->dispIn[n-1]);
        int yIn1  = mapY(d->dispIn[n]);
        int yOut0 = mapY(d->dispOut[n-1]);
        int yOut1 = mapY(d->dispOut[n]);

        NT_drawShapeI(kNT_line, x0, yIn0,  x1, yIn1,  6);
        NT_drawShapeI(kNT_line, x0, yOut0, x1, yOut1, 15);
    }

    int lockCol = d->locked ? 15 : 4;
    NT_drawShapeI(kNT_rectangle, width-8, 2, width-2, 8, lockCol);

    return true;
}

static uint32_t hasCustomUi(_NT_algorithm *) { return 0; }
static void customUi(_NT_algorithm *, const _NT_uiData &) {}

static void calculateStaticRequirements(_NT_staticRequirements& req) { req.dram = 0; }
static void initialise(_NT_staticMemoryPtrs& ptrs, const _NT_staticRequirements& req) {
    (void)ptrs; (void)req;
}

static void calculateRequirements(_NT_algorithmRequirements &req, const int32_t* specifications) {
    (void)specifications;
    req.numParameters = kNumParams;
    req.sram = sizeof(_CessnaAwg);
    req.dram = 0;
    req.dtc  = sizeof(_CessnaAwg_DTC);
    req.itc  = 0;
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