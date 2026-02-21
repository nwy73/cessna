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
 * - Outputs: OUT1 = resynth audio, OUT2 = error (input - resynth).
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

    // One-cycle accumulation for bin means
    float binSum[kBins];
    int   binCount[kBins];
    int   cycleSampleCount;
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
static void calculateRequirements(_NT_algorithmRequirements &req) {
    req.sram = sizeof(_CessnaAwg);
    req.dtc  = sizeof(_CessnaAwg_DTC);
}

static _NT_algorithm* construct(const _NT_algorithmDescriptor*, void *sram, void *dtc) {
    auto *d = new(dtc) _CessnaAwg_DTC();
    auto *self = new(sram) _CessnaAwg(d);

    d->sampleRate = (NT_globals.sampleRate > 0) ? (float)NT_globals.sampleRate : kDefaultSampleRate;

    d->depthD = 0.5f;
    d->adaptAlpha = 0.2f;
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

    // Clear accumulators
    for (int i=0;i<kBins;i++){ d->binSum[i]=0.0f; d->binCount[i]=0; }
    d->cycleSampleCount = 0;

    d->phase = 0.0f;
    d->phaseInc = 1.0f / d->smoothPeriod;

    // Display init
    for (int i=0;i<kDisplayPoints;i++){ d->dispIn[i]=0.0f; d->dispOut[i]=0.0f; }
    d->dispInit = true;

    return self;
}

// ----------------------------
// Parameter changes
// ----------------------------
static void parameterChanged(_NT_algorithm *base, int param, int value) {
    auto *self = (_CessnaAwg*)base;
    if (!self || !self->dtc) return;
    auto *d = self->dtc;

    // NOTE: The NT parameter value encoding varies (int16, etc.). This example assumes
    // parameters arrive as int with useful scaling handled by the descriptor.
    // In real plug-ins you typically convert based on parameter type/range.

    if (param == kParamDepth) {
        d->depthD = (float)value / 1000.0f; // placeholder scaling
    } else if (param == kParamAdapt) {
        d->adaptAlpha = (float)value / 1000.0f;
    } else if (param == kParamCurve) {
        d->curve = (float)value / 1000.0f;
    } else if (param == kParamLock) {
        d->lockStrength = (float)value / 1000.0f;
    } else if (param == kParamFreeze) {
        d->freeze = (value != 0);
    } else if (param == kParamCapture) {
        if (value != 0) d->pendingCapture = true; // momentary trigger
    } else if (param >= kParamX1 && param <= kParamX16) {
        int i = param - kParamX1;
        // Expect symmetric mapping already produces bipolar range; placeholder scaling:
        d->X[i] = (float)value / 1000.0f;
    }
}

// ----------------------------
// Core: per-sample resynth evaluation
// ----------------------------
static inline float evalLinear(const float *p, float phase01) {
    float x = phase01 * (float)kBins;
    int i = (int)floorf(x);
    float t = x - (float)i;
    i &= 15;
    int i1 = (i+1) & 15;
    return lerpf(p[i], p[i1], t);
}

// ----------------------------
// Cycle detector + bin accumulation
// ----------------------------
static inline void resetAccumulators(_CessnaAwg_DTC *d) {
    for (int i=0;i<kBins;i++){ d->binSum[i]=0.0f; d->binCount[i]=0; }
    d->cycleSampleCount = 0;
}

static inline void finalizeCycle(_CessnaAwg_DTC *d) {
    // Compute new bin means
    for (int i=0;i<kBins;i++) {
        if (d->binCount[i] > 0) d->newA[i] = d->binSum[i] / (float)d->binCount[i];
        else d->newA[i] = d->A[i];
    }

    // Capture logic: only if locked, and capture armed
    if (d->pendingCapture && d->locked) {
        for (int i=0;i<kBins;i++) d->A[i] = d->newA[i];
        d->freeze = true;          // Capture => Freeze ON
        d->pendingCapture = false; // disarm
    } else {
        // Normal update (unless frozen)
        float alpha = d->freeze ? 0.0f : clampf(d->adaptAlpha, 0.0f, 1.0f);
        for (int i=0;i<kBins;i++) {
            d->A[i] += alpha * (d->newA[i] - d->A[i]);
        }
    }

    // Update S[i] with deviations (pre-interpolation)
    float D = d->depthD;
    for (int i=0;i<kBins;i++) {
        d->S[i] = d->A[i] + D * d->X[i];
    }

    // Recompute curved slopes
    d->herm.computeSlopes(d->S);

    // Update phase increment based on smoothed period
    if (d->smoothPeriod > 4.0f) {
        d->phaseInc = 1.0f / d->smoothPeriod;
    }

    // Update display buffers by resampling one cycle to kDisplayPoints.
    // We use the most recent cycle's input bins as a crude display of input;
    // in practice you'd resample the raw captured cycle samples.
    // For v1 we approximate input display using newA (bin means) interpolated.
    for (int n=0;n<kDisplayPoints;n++) {
        float ph = (float)n / (float)kDisplayPoints;
        // Input display from newA (most recent analysis)
        float inLin = evalLinear(d->newA, ph);
        float inCur = d->herm.eval(d->newA, ph); // slopes not computed for newA; harmless but not ideal
        (void)inCur;
        d->dispIn[n] = inLin;

        // Output display from S
        float outLin = evalLinear(d->S, ph);
        float outCur = d->herm.eval(d->S, ph);
        d->dispOut[n] = lerpf(outLin, outCur, d->curve);
    }

    // Reset for next cycle
    resetAccumulators(d);
}

static inline bool processCycleDetection(_CessnaAwg_DTC *d, float x) {
    // Returns true if a cycle boundary was detected.
    float xd = d->dc.process(x);

    // Envelope follower for threshold
    float ax = fabsf(xd);
    d->env += kEnvBeta * (ax - d->env);
    d->thresh = fmaxf(kThreshMin, kThreshK * d->env);

    float T = d->thresh;

    // Track state
    if (xd < -T) d->wasBelowNeg = true;

    bool edge = false;
    if (d->wasBelowNeg && xd > +T) {
        // Candidate rising crossing
        int period = d->samplesSinceEdge;
        if (period >= d->minPeriodSamps && period <= d->maxPeriodSamps) {
            d->lastPeriodSamples = period;
            // Smooth period
            float gamma = 0.15f;
            d->smoothPeriod += gamma * ((float)period - d->smoothPeriod);
            d->locked = true;
            edge = true;
        } else {
            // Too fast -> ignore, too slow -> may be unlock
            if (period > d->maxPeriodSamps) d->locked = false;
        }
        d->samplesSinceEdge = 0;
        d->wasBelowNeg = false;
    } else {
        d->samplesSinceEdge++;
        if (d->samplesSinceEdge > d->maxPeriodSamps * 2) {
            d->locked = false;
            // In unlocked state, do not auto-clear wasBelowNeg; keep it conservative.
            d->samplesSinceEdge = d->maxPeriodSamps * 2;
        }
    }

    return edge;
}

static inline void accumulateBin(_CessnaAwg_DTC *d, float x) {
    // Assign sample within current cycle to a bin by proportional position.
    // We don't know final period until edge; we approximate using smoothPeriod.
    // This works well when locked; when not locked we still accumulate but will likely freeze.
    int P = (int)fmaxf(4.0f, d->smoothPeriod);
    int idx = d->cycleSampleCount;
    if (idx < 0) idx = 0;
    if (idx >= P) idx = P-1;

    int bin = (idx * kBins) / P; // 0..15
    if (bin < 0) bin = 0;
    if (bin > 15) bin = 15;

    d->binSum[bin] += x;
    d->binCount[bin] += 1;
    d->cycleSampleCount++;
}

// ----------------------------
// step – audio processing
// ----------------------------
static void step(_NT_algorithm *base, float **inBuf, float **outBuf, int numFrames) {
    auto *self = (_CessnaAwg*)base;
    if (!self || !self->dtc) return;
    auto *d = self->dtc;

    float *in  = (inBuf && inBuf[0]) ? inBuf[0] : nullptr;
    float *out1 = (outBuf && outBuf[0]) ? outBuf[0] : nullptr;
    float *out2 = (outBuf && outBuf[1]) ? outBuf[1] : nullptr;

    if (!in) {
        // No input: just output resynth based on current bins
        for (int n=0;n<numFrames;n++) {
            // Build S on the fly
            float D = d->depthD;
            for (int i=0;i<kBins;i++) d->S[i] = d->A[i] + D*d->X[i];
            d->herm.computeSlopes(d->S);

            float lin = evalLinear(d->S, d->phase);
            float cur = d->herm.eval(d->S, d->phase);
            float y = lerpf(lin, cur, d->curve);

            if (out1) out1[n] = y;
            if (out2) out2[n] = 0.0f;

            d->phase += d->phaseInc;
            if (d->phase >= 1.0f) d->phase -= 1.0f;
        }
        return;
    }

    for (int n=0;n<numFrames;n++) {
        float x = in[n];

        // Tick highlight on change of X[i] (control-rate-ish, cheap)
        for (int i=0;i<kBins;i++) {
            float dx = fabsf(d->X[i] - d->X_prev[i]);
            if (dx > kTickChangeThresh) d->tickHi[i] = 1.0f;
            d->X_prev[i] = d->X[i];
        }

        // Cycle detection
        bool edge = processCycleDetection(d, x);

        // Accumulate input samples into bins for the current cycle
        accumulateBin(d, x);

        if (edge) {
            // If lock strength == 1, hard reset phase at cycle boundary
            if (d->lockStrength >= 0.999f) {
                d->phase = 0.0f;
            }

            finalizeCycle(d);
        }

        // Compute resynth sample from S (updated once per cycle, but we ensure S exists)
        // Build S lazily if not yet built
        if (d->cycleSampleCount == 0) {
            float D = d->depthD;
            for (int i=0;i<kBins;i++) d->S[i] = d->A[i] + D*d->X[i];
            d->herm.computeSlopes(d->S);
        }

        float yLin = evalLinear(d->S, d->phase);
        float yCur = d->herm.eval(d->S, d->phase);
        float y = lerpf(yLin, yCur, d->curve);

        // Outputs
        if (out1) out1[n] = y;
        if (out2) out2[n] = x - y; // raw error

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
    NT_drawShapeI(kNT_rect, width-8, 2, width-2, 8, lockCol);

    return true;
}

static uint32_t hasCustomUi(_NT_algorithm *) {
    return 0; // no custom controls (yet)
}

static void customUi(_NT_algorithm *, const _NT_uiData &) {}

// ----------------------------
// Descriptor
// ----------------------------
static _NT_algorithmDescriptor descriptor = {
    .apiVersion = kNT_apiVersion,
    .name = "Cessna++ AWG",
    .flags = 0,
    .calculateRequirements = calculateRequirements,
    .construct = construct,
    .destruct = nullptr,
    .parameterChanged = parameterChanged,
    .step = step,
    .draw = draw,
    .hasCustomUi = hasCustomUi,
    .customUi = customUi,
};

// ----------------------------
// Plug-in entry point
// ----------------------------
extern "C" {
_NT_pluginInfo NT_pluginInfo = {
    .apiVersion = kNT_apiVersion,
    .pluginName = "Cessna++",
    .numAlgorithms = 1,
    .algorithms = &descriptor,
};
}
