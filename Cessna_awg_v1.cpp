/*
 * Cessna AWG (v1) — CV-driven Arbitrary Waveform Generator
 * ---------------------------------------------------------
 * A faithful recreation of James Cessna's 1970 AWG design:
 * a wavetable oscillator driven by a 1V/oct pitch CV input,
 * with 16 manually programmable amplitude bins defining the
 * waveform shape. 0V = C4 (261.63 Hz).
 *
 * Extensions beyond the original:
 * - Smoothing: linear <-> cubic Hermite interpolation blend
 * - Modulation: Ripple, Sequencer, and Drunk modes for
 *   animating the bin values in real time
 * - Depth: scales overall output amplitude
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
// Constants
// ----------------------------
static constexpr int   kBins              = 16;
static constexpr int   kDisplayPoints     = 128;
static constexpr float kDefaultSampleRate = 48000.0f;
static constexpr float kC4Hz             = 261.63f;   // 0V = C4

// ----------------------------
// Helpers
// ----------------------------
static inline float clampf(float x, float lo, float hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}
static inline float lerpf(float a, float b, float t) {
    return a + (b - a) * t;
}

// LCG random number generator (same as v39)
static inline float lcgRand(uint32_t &state) {
    state = state * 1664525u + 1013904223u;
    return ((float)(state >> 1) / (float)0x7fffffff) - 1.0f; // -1..1
}

// ----------------------------
// Forward declarations
// ----------------------------
static inline float evalLinear(const float *bins, float phase01);
static void seqBuildPattern(struct _CessnaAwg_DTC *d);
static void computeModX(struct _CessnaAwg_DTC *d);

// ----------------------------
// Cubic Hermite interpolator
// ----------------------------
struct PeriodicHermite16 {
    float m[kBins] = {0};

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
            if (dp == 0.0f) { m[i] = 0.0f; m[(i+1)&15] = 0.0f; continue; }
            if ((m[i] / dp) < 0.0f)        m[i]        = 0.0f;
            if ((m[(i+1)&15] / dp) < 0.0f) m[(i+1)&15] = 0.0f;
        }
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
        float p0 = p[i], p1 = p[i1];
        float m0 = m[i], m1 = m[i1];
        float t2 = t*t, t3 = t2*t;
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

    // Waveform bins — direct amplitude values, ±5V
    float B[kBins];     // base bin values (from user parameters)
    float S[kBins];     // combined: B[i] + modX[i], used for output

    float depthD;       // 0..1 output amplitude scale
    float curve;        // 0..1 linear→hermite blend
    bool  stepped;      // true = no interpolation, hold each bin value

    float phase;
    float phaseInc;
    float phaseInc;

    PeriodicHermite16 herm;

    // --- Modulation state ---
    int   modMode;          // 0=Off 1=Ripple 2=Sequencer 3=Drunk
    float modAmt;
    float rippleRate;
    float lfoPhase;
    float lfoPhaseInc;

    int   seqSpeed;
    int   modSampleCount;   // counts samples between seq/drunk advances
    int   modSamplePeriod;  // samples per advance step (derived from Seq Rate)
    float seqSmooth;
    float seqDrift;
    int   seqPattern;
    float seqTarget[kBins];
    float seqCurrent[kBins];
    int   seqCycleCount;
    uint32_t seqRng;

    float drunkWander;
    float drunkValues[kBins];
    uint32_t drunkRng;

    float modX[kBins];      // modulation contribution, ±1 (scaled by modAmt)

    // Waveform slots for morph
    float slotA[kBins];     // stored waveform A (locked when saveA=true)
    float slotB[kBins];     // stored waveform B (locked when saveB=true)
    bool  saveA;            // true = slotA is locked; false = slotA tracks B[]
    bool  saveB;            // true = slotB is locked; false = slotB tracks B[]
    float morph;            // 0..1 crossfade between slotA and slotB

    // Display
    float dispWave[kDisplayPoints];

    // Grayout / UI
    bool  needsGrayOutUpdate;
    int   algorithmIdx;

    _CessnaAwg_DTC() {
        sampleRate = kDefaultSampleRate;

        depthD  = 1.0f;
        curve   = 0.0f;
        stepped = false;

        phase          = 0.0f;
        phaseInc       = kC4Hz / kDefaultSampleRate;

        modMode     = 0;
        modAmt      = 0.5f;
        rippleRate  = 0.1f;
        lfoPhase    = 0.0f;
        lfoPhaseInc = rippleRate / kDefaultSampleRate;

        seqSpeed        = 4;
        modSampleCount  = 0;
        // Default Seq Rate param=29 → seqSpeed=4 → ~83ms at 48kHz
        modSamplePeriod = (int)(kDefaultSampleRate * 0.083f);
        seqSmooth     = 0.5f;
        seqDrift      = 0.2f;
        seqPattern    = 0;
        seqCycleCount = 0;
        seqRng        = 0x12345678u;

        drunkWander = 0.05f;
        drunkRng    = 0x87654321u;

        for (int i=0;i<kBins;i++) {
            B[i]           = 0.0f;   // silence default
            S[i]           = 0.0f;
            seqTarget[i]   = 0.0f;
            seqCurrent[i]  = 0.0f;
            drunkValues[i] = 0.0f;
            modX[i]        = 0.0f;
        }

        // Initialise sequencer pattern so seq mode has something to play immediately
        // (seqBuildPattern needs seqPattern, seqDrift, seqRng already set above)
        {
            // Inline ramp init matching case 0
            for (int i=0;i<kBins;i++)
                seqTarget[i] = ((float)i / (float)(kBins-1)) * 2.0f - 1.0f;
            // Copy to current so there's no initial slew from zero
            for (int i=0;i<kBins;i++)
                seqCurrent[i] = seqTarget[i];
        }

        // Seed drunk values with a small random spread so it's immediately active
        {
            uint32_t rng = drunkRng;
            for (int i=0;i<kBins;i++) {
                rng = rng * 1664525u + 1013904223u;
                drunkValues[i] = ((float)(rng >> 1) / (float)0x7fffffff - 1.0f) * 0.3f;
            }
            drunkRng = rng;
        }

        for (int i=0;i<kDisplayPoints;i++) dispWave[i] = 0.0f;

        saveA = false;
        saveB = false;
        morph = 0.0f;
        for (int i=0;i<kBins;i++) { slotA[i] = 0.0f; slotB[i] = 0.0f; }

        needsGrayOutUpdate = true;
        algorithmIdx = -1;
    }
};

// ----------------------------
// Plugin struct
// ----------------------------
struct _CessnaAwg : public _NT_algorithm {
    _CessnaAwg_DTC *dtc;
    _CessnaAwg(_CessnaAwg_DTC *d) : dtc(d) {}
};

// ----------------------------
// Parameter enum
// ----------------------------
enum {
    kParamCVIn = 0,
    kParamOut1, kParamOut1Mode,

    kParamDepth,
    kParamCurve,
    kParamStepped,

    // Bins (waveform page)
    kParamB1,
    kParamB16 = kParamB1 + 15,

    // Mod page
    kParamModMode,
    kParamModAmt,
    kParamRippleRate,
    kParamSeqSpeed,
    kParamSeqSmooth,
    kParamSeqDrift,
    kParamSeqPattern,
    kParamDrunkWander,

    // Morph
    kParamSaveA,
    kParamSaveB,
    kParamMorph,
};
static constexpr int kNumParams = kParamMorph + 1;

// ----------------------------
// String tables
// ----------------------------
static const char* modModeStrings[] = {"Off", "Ripple", "Sequencer", "Drunk", nullptr};
static const char* seqPatStrings[]  = {"Ramp", "Triangle", "Random", "Alternating", nullptr};
static const char* onOffStrings[]   = {"Off", "On", nullptr};

// ----------------------------
// Parameter descriptors
// ----------------------------
static const _NT_parameter g_parameters[kNumParams] = {
    NT_PARAMETER_CV_INPUT("CV in", 0, 1)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("out", 0, 13)

    { .name="Level",    .min=0,  .max=100, .def=100, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Smoothing",.min=0,  .max=100, .def=0,   .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Stepped",  .min=0,  .max=1,   .def=0,   .unit=kNT_unitEnum,    .scaling=kNT_scalingNone, .enumStrings=onOffStrings },

    // Bins 1–16: amplitude ±100% = ±5V
    { .name="Bin 1",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 2",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 3",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 4",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 5",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 6",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 7",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 8",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 9",  .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 10", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 11", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 12", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 13", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 14", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 15", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 16", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },

    // Mod page
    { .name="Mod Mode",    .min=0,.max=3,   .def=0,  .unit=kNT_unitEnum,    .scaling=kNT_scalingNone, .enumStrings=modModeStrings },
    { .name="Mod Amt",     .min=0,.max=100, .def=50, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Ripple Rate", .min=1,.max=1000,.def=10, .unit=kNT_unitNone,    .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Seq Rate",    .min=1,.max=32,  .def=29, .unit=kNT_unitNone,    .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Seq Smooth",  .min=0,.max=100, .def=50, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Seq Drift",   .min=0,.max=100, .def=20, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
    { .name="Seq Pattern", .min=0,.max=3,   .def=0,  .unit=kNT_unitEnum,    .scaling=kNT_scalingNone, .enumStrings=seqPatStrings },
    { .name="Wander",      .min=0,.max=100, .def=20, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },

    // Morph page
    { .name="Save A", .min=0,.max=1,.def=0, .unit=kNT_unitEnum, .scaling=kNT_scalingNone, .enumStrings=onOffStrings },
    { .name="Save B", .min=0,.max=1,.def=0, .unit=kNT_unitEnum, .scaling=kNT_scalingNone, .enumStrings=onOffStrings },
    { .name="Morph",  .min=0,.max=100,.def=0, .unit=kNT_unitPercent, .scaling=kNT_scalingNone, .enumStrings=nullptr },
};

// ----------------------------
// Page index arrays
// ----------------------------
static const uint8_t g_pageRoutingIdx[] = {
    (uint8_t)kParamCVIn,
    (uint8_t)kParamOut1, (uint8_t)kParamOut1Mode,
};

static const uint8_t g_pageMainIdx[] = {
    (uint8_t)kParamDepth,
    (uint8_t)kParamCurve,
    (uint8_t)kParamStepped,
};

static const uint8_t g_pageBinsIdx[] = {
    (uint8_t)kParamB1,      (uint8_t)(kParamB1+1),  (uint8_t)(kParamB1+2),  (uint8_t)(kParamB1+3),
    (uint8_t)(kParamB1+4),  (uint8_t)(kParamB1+5),  (uint8_t)(kParamB1+6),  (uint8_t)(kParamB1+7),
    (uint8_t)(kParamB1+8),  (uint8_t)(kParamB1+9),  (uint8_t)(kParamB1+10), (uint8_t)(kParamB1+11),
    (uint8_t)(kParamB1+12), (uint8_t)(kParamB1+13), (uint8_t)(kParamB1+14), (uint8_t)(kParamB1+15),
};

static const uint8_t g_pageModIdx[] = {
    (uint8_t)kParamModMode,
    (uint8_t)kParamModAmt,
    (uint8_t)kParamRippleRate,
    (uint8_t)kParamSeqSpeed,
    (uint8_t)kParamSeqSmooth,
    (uint8_t)kParamSeqDrift,
    (uint8_t)kParamSeqPattern,
    (uint8_t)kParamDrunkWander,
};

static const uint8_t g_pageMorphIdx[] = {
    (uint8_t)kParamSaveA,
    (uint8_t)kParamSaveB,
    (uint8_t)kParamMorph,
};

static const _NT_parameterPage g_pages_arr[] = {
    { .name="Routing", .numParams=(uint8_t)ARRAY_SIZE(g_pageRoutingIdx), .group=0, .unused={0,0}, .params=g_pageRoutingIdx },
    { .name="Main",    .numParams=(uint8_t)ARRAY_SIZE(g_pageMainIdx),    .group=0, .unused={0,0}, .params=g_pageMainIdx },
    { .name="Waveform",.numParams=(uint8_t)ARRAY_SIZE(g_pageBinsIdx),    .group=1, .unused={0,0}, .params=g_pageBinsIdx },
    { .name="Mod",     .numParams=(uint8_t)ARRAY_SIZE(g_pageModIdx),     .group=0, .unused={0,0}, .params=g_pageModIdx },
    { .name="Morph",   .numParams=(uint8_t)ARRAY_SIZE(g_pageMorphIdx),   .group=0, .unused={0,0}, .params=g_pageMorphIdx },
};

static const _NT_parameterPages g_pages = {
    .numPages = ARRAY_SIZE(g_pages_arr),
    .pages    = g_pages_arr,
};

// ----------------------------
// Helper: get parameter as 0..1
// ----------------------------
#define getPct01(p)      ((float)self->v[(p)] * 0.01f)
#define getPctSigned(p)  ((float)self->v[(p)] * 0.01f * 5.0f)  // ±100% → ±5V

// ----------------------------
// Linear wavetable evaluation
// ----------------------------
static inline float evalLinear(const float *bins, float phase01) {
    float x  = phase01 * (float)kBins;
    int   i  = (int)floorf(x) & 15;
    float t  = x - floorf(x);
    return lerpf(bins[i], bins[(i+1)&15], t);
}

static inline float evalStepped(const float *bins, float phase01) {
    int i = (int)floorf(phase01 * (float)kBins) & 15;
    return bins[i];
}

// ----------------------------
// combineS: merge B[] + modX[] into S[]
// ----------------------------
static void combineS(struct _CessnaAwg_DTC *d) {
    // Update slots from live bins if not saved
    if (!d->saveA) for (int i=0;i<kBins;i++) d->slotA[i] = d->B[i];
    if (!d->saveB) for (int i=0;i<kBins;i++) d->slotB[i] = d->B[i];

    float D = d->depthD;
    for (int i=0;i<kBins;i++) {
        // Morph between the two slots
        float base = lerpf(d->slotA[i], d->slotB[i], d->morph);
        float v = D * clampf(base + d->modX[i] * 5.0f, -5.0f, 5.0f);
        d->S[i] = clampf(v, -5.0f, 5.0f);
    }
}

// ----------------------------
// Modulation
// ----------------------------
static void seqBuildPattern(struct _CessnaAwg_DTC *d) {
    switch (d->seqPattern) {
        case 0: // Ramp
            for (int i=0;i<kBins;i++)
                d->seqTarget[i] = ((float)i / (float)(kBins-1)) * 2.0f - 1.0f;
            break;
        case 1: { // Triangle
            for (int i=0;i<kBins;i++) {
                float t = (float)i / (float)kBins;
                d->seqTarget[i] = (t < 0.5f) ? (t*4.0f - 1.0f) : (3.0f - t*4.0f);
            }
            break;
        }
        case 2: // Random
            for (int i=0;i<kBins;i++)
                d->seqTarget[i] = clampf(d->seqTarget[i] + lcgRand(d->seqRng) * d->seqDrift, -1.0f, 1.0f);
            break;
        case 3: // Alternating
            for (int i=0;i<kBins;i++)
                d->seqTarget[i] = (i & 1) ? -1.0f : 1.0f;
            break;
    }
    if (d->seqPattern != 2 && d->seqDrift > 0.0f) {
        for (int i=0;i<kBins;i++)
            d->seqTarget[i] = clampf(d->seqTarget[i] + lcgRand(d->seqRng) * d->seqDrift, -1.0f, 1.0f);
    }
}

static void seqRotate(struct _CessnaAwg_DTC *d) {
    float tmp = d->seqTarget[kBins-1];
    for (int i=kBins-1;i>0;i--) d->seqTarget[i] = d->seqTarget[i-1];
    d->seqTarget[0] = tmp;
}

static void seqAdvance(struct _CessnaAwg_DTC *d) {
    if (d->seqPattern == 2) {
        seqBuildPattern(d);
    } else {
        seqRotate(d);
        if (d->seqDrift > 0.0f) {
            for (int i=0;i<kBins;i++)
                d->seqTarget[i] = clampf(d->seqTarget[i] + lcgRand(d->seqRng) * d->seqDrift, -1.0f, 1.0f);
        }
    }
}

static void drunkAdvance(struct _CessnaAwg_DTC *d) {
    if (d->drunkWander < 0.001f) {
        for (int i=0;i<kBins;i++) d->drunkValues[i] *= 0.9f;
        return;
    }
    for (int i=0;i<kBins;i++) {
        float step = lcgRand(d->drunkRng) * d->drunkWander;
        d->drunkValues[i] = clampf(d->drunkValues[i] + step, -1.0f, 1.0f);
    }
}

static void computeModX(struct _CessnaAwg_DTC *d) {
    switch (d->modMode) {
        case 0: // Off
            for (int i=0;i<kBins;i++) d->modX[i] = 0.0f;
            break;
        case 1: { // Ripple
            d->lfoPhase += d->lfoPhaseInc;
            if (d->lfoPhase >= 1.0f) d->lfoPhase -= 1.0f;
            for (int i=0;i<kBins;i++) {
                float binPhase = d->lfoPhase + (float)i / (float)kBins;
                if (binPhase >= 1.0f) binPhase -= 1.0f;
                d->modX[i] = d->modAmt * sinf(2.0f * M_PI_F * binPhase);
            }
            break;
        }
        case 2: { // Sequencer
            float slew = 1.0f - d->seqSmooth * 0.999f;
            for (int i=0;i<kBins;i++) {
                d->seqCurrent[i] += slew * (d->seqTarget[i] - d->seqCurrent[i]);
                d->modX[i] = d->modAmt * d->seqCurrent[i];
            }
            break;
        }
        case 3: // Drunk
            for (int i=0;i<kBins;i++)
                d->modX[i] = d->modAmt * d->drunkValues[i];
            break;
    }
}

// ----------------------------
// Grayout
// ----------------------------
static void updateModGrayOut(_NT_algorithm* base, int modeOverride = -1) {
    auto *self = (_CessnaAwg*)base;
    
    auto *d = self->dtc;
    if (!d) return;

    if (d->algorithmIdx < 0) {
        d->needsGrayOutUpdate = true;
        return;
    }
    uint32_t idx = (uint32_t)d->algorithmIdx;

    int mode = (modeOverride >= 0) ? modeOverride : d->modMode;

    bool isOff  = (mode == 0);
    bool ripple = (mode == 1);
    bool seq    = (mode == 2);
    bool drunk  = (mode == 3);
    bool stepped = (self->v[kParamStepped] != 0);

    NT_setParameterGrayedOut(idx, kParamCurve,       stepped);
    NT_setParameterGrayedOut(idx, kParamModAmt,      isOff);
    NT_setParameterGrayedOut(idx, kParamRippleRate,  !ripple);
    NT_setParameterGrayedOut(idx, kParamSeqSpeed,    !seq);
    NT_setParameterGrayedOut(idx, kParamSeqSmooth,   !seq);
    NT_setParameterGrayedOut(idx, kParamSeqDrift,    !seq);
    NT_setParameterGrayedOut(idx, kParamSeqPattern,  !seq);
    NT_setParameterGrayedOut(idx, kParamDrunkWander, !drunk);
}

// ----------------------------
// parameterChanged()
// ----------------------------
static void parameterChanged(_NT_algorithm *base, int p) {
    auto *self = (_CessnaAwg*)base;
    
    auto *d = self->dtc;
    if (!d) return;

    switch (p) {
        case kParamDepth:   d->depthD  = getPct01(kParamDepth); break;  // Level
        case kParamCurve:   d->curve   = getPct01(kParamCurve); break;
        case kParamStepped: d->stepped = (self->v[kParamStepped] != 0); break;

        case kParamModMode: d->modMode = self->v[kParamModMode]; break;
        case kParamModAmt:  d->modAmt  = getPct01(kParamModAmt); break;
        case kParamRippleRate: {
            d->rippleRate  = (float)self->v[kParamRippleRate] * 0.01f;
            d->lfoPhaseInc = d->rippleRate / d->sampleRate;
            break;
        }
        case kParamSeqSpeed: {
            // param 1..32, higher = faster
            // Map to sample period: param=1 → ~2s, param=32 → ~50ms
            float pct = (float)self->v[kParamSeqSpeed] / 32.0f; // 0..1
            // Exponential mapping for musical feel
            float periodSecs = 2.0f * powf(0.025f / 2.0f, pct); // 2s → 0.025s
            d->modSamplePeriod = (int)(d->sampleRate * periodSecs);
            if (d->modSamplePeriod < 1) d->modSamplePeriod = 1;
            break;
        }
        case kParamSeqSmooth:   d->seqSmooth  = getPct01(kParamSeqSmooth); break;
        case kParamSeqDrift:    d->seqDrift   = getPct01(kParamSeqDrift); break;
        case kParamSeqPattern:  d->seqPattern = self->v[kParamSeqPattern]; seqBuildPattern(d); break;
        case kParamDrunkWander: d->drunkWander = getPct01(kParamDrunkWander); break;

        case kParamSaveA: {
            d->saveA = (self->v[kParamSaveA] != 0);
            if (d->saveA) {
                // Lock slot A to current bin values
                for (int i=0;i<kBins;i++) d->slotA[i] = d->B[i];
            }
            break;
        }
        case kParamSaveB: {
            d->saveB = (self->v[kParamSaveB] != 0);
            if (d->saveB) {
                // Lock slot B to current bin values
                for (int i=0;i<kBins;i++) d->slotB[i] = d->B[i];
            }
            break;
        }
        case kParamMorph: d->morph = getPct01(kParamMorph); break;

        default:
            if (p >= kParamB1 && p <= kParamB16) {
                int i = p - kParamB1;
                d->B[i] = getPctSigned(p);  // ±100% → ±5V
            }
            break;
    }

    // Always defer grayout to step() — self->v updated by NT after callback
    d->needsGrayOutUpdate = true;
}

// ----------------------------
// updateDisplay: draw waveform from S[]
// ----------------------------
static void updateDisplay(_CessnaAwg_DTC *d) {
    for (int k=0;k<kDisplayPoints;k++) {
        float ph = (float)k / (float)(kDisplayPoints - 1);
        float y;
        if (d->stepped) {
            y = evalStepped(d->S, ph);
        } else {
            float lin = evalLinear(d->S, ph);
            float cur = (d->curve > 0.0001f) ? d->herm.eval(d->S, ph) : lin;
            y = lerpf(lin, cur, d->curve);
        }
        d->dispWave[k] = y;
    }
}

// ----------------------------
// step()
// ----------------------------
static void step(_NT_algorithm *base, float *busFrames, int numFramesBy4) {
    auto *self = (_CessnaAwg*)base;
    auto *d = self->dtc;
    if (!d) return;

    const int numFrames = numFramesBy4 * 4;

    // Cache algorithm index on first step
    if (d->algorithmIdx < 0) {
        d->algorithmIdx = (int)NT_algorithmIndex(self);
    }

    // Deferred grayout
    if (d->needsGrayOutUpdate) {
        d->needsGrayOutUpdate = false;
        d->modMode  = self->v[kParamModMode];
        updateModGrayOut(base);
    }

    const int cvBus  = self->v[kParamCVIn];
    const int outBus = self->v[kParamOut1];
    const int outMode= self->v[kParamOut1Mode];

    // busPtr helper
    auto busPtr = [&](int busVal) -> float* {
        if (!busFrames || busVal <= 0) return nullptr;
        int busIndex = busVal - 1;
        if (busIndex < 0 || busIndex >= 28) return nullptr;
        return busFrames + busIndex * numFrames;
    };

    float *cv  = busPtr(cvBus);
    float *out = busPtr(outBus);

    for (int n=0;n<numFrames;n++) {
        // --- Pitch from CV ---
        float cvVal = cv ? cv[n] : 0.0f;
        float freq  = kC4Hz * powf(2.0f, cvVal);
        d->phaseInc = freq / d->sampleRate;

        // --- Modulation ---
        computeModX(d);

        // --- Combine bins + mod ---
        combineS(d);
        if (d->curve > 0.0001f) d->herm.computeSlopes(d->S);

        // --- Evaluate wavetable ---
        float y;
        if (d->stepped) {
            y = evalStepped(d->S, d->phase);
        } else {
            float lin = evalLinear(d->S, d->phase);
            float cur = (d->curve > 0.0001f) ? d->herm.eval(d->S, d->phase) : lin;
            y = lerpf(lin, cur, d->curve);
        }

        // --- Output ---
        if (out) {
            if (outMode == 0) out[n] += y;
            else              out[n]  = y;
        }

        // --- Advance phase ---
        d->phase += d->phaseInc;
        if (d->phase >= 1.0f) d->phase -= 1.0f;

        // --- Seq/Drunk advance on sample timer (pitch-independent) ---
        if (d->modMode == 2 || d->modMode == 3) {
            d->modSampleCount++;
            if (d->modSampleCount >= d->modSamplePeriod) {
                d->modSampleCount = 0;
                if (d->modMode == 2) seqAdvance(d);
                if (d->modMode == 3) drunkAdvance(d);
            }
        }
    }

    // Update display waveform
    updateDisplay(d);
}

// ----------------------------
// draw()
// ----------------------------
static bool draw(_NT_algorithm *base) {
    auto *self = (_CessnaAwg*)base;
    auto *d = self->dtc;
    if (!d) return false;

    const int W = 256;
    const int H = 64;

    int top    = 11;
    int bottom = H - 2;
    int plotH  = bottom - top;

    auto mapY = [&](float v)->int {
        float vv = clampf(v, -5.0f, 5.0f);
        float yn = 0.5f * (1.0f - (vv / 5.0f));
        int y = top + (int)roundf(yn * (float)plotH);
        if (y < top)  y = top;
        if (y >= H)   y = H-1;
        return y;
    };

    // Draw waveform from dispWave[]
    for (int n=1;n<kDisplayPoints;n++) {
        int x0 = (int)roundf(((float)(n-1) / (float)(kDisplayPoints-1)) * (float)(W-1));
        int x1 = (int)roundf(((float)n     / (float)(kDisplayPoints-1)) * (float)(W-1));
        int y0 = mapY(d->dispWave[n-1]);
        int y1 = mapY(d->dispWave[n]);
        NT_drawShapeI(kNT_line, x0, y0, x1, y1, 15);
        NT_drawShapeI(kNT_line, x0, y0, x1, y1, 15); // double pass = bold
    }

    // Draw bin tick marks at the bottom
    for (int i=0;i<kBins;i++) {
        float ph = ((float)i + 0.5f) / (float)kBins;
        int x = (int)roundf(ph * (float)(W-1));
        NT_drawShapeI(kNT_line, x, H-4, x, H-1, 6);
    }

    return false; // let NT draw the standard parameter line at top
}

// ----------------------------
// NT API registration
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
    (void)req; (void)specifications;

    auto *d    = new(ptrs.dtc)  _CessnaAwg_DTC();
    auto *self = new(ptrs.sram) _CessnaAwg(d);

    self->parameters    = g_parameters;
    self->parameterPages = &g_pages;

    float sr = (float)NT_globals.sampleRate;
    if (sr < 8000.0f || sr > 192000.0f) sr = kDefaultSampleRate;
    d->sampleRate  = sr;
    d->lfoPhaseInc = d->rippleRate / sr;
    d->modSamplePeriod = (int)(sr * 0.083f);
    d->phaseInc    = kC4Hz / sr;

    return self;
}

static uint32_t hasCustomUi(_NT_algorithm *) { return 0; }
static void customUi(_NT_algorithm *, const _NT_uiData &) {}
static void calculateStaticRequirements(_NT_staticRequirements& req) { req.dram = 0; }
static void initialise(_NT_staticMemoryPtrs& ptrs, const _NT_staticRequirements& req) {
    (void)ptrs; (void)req;
}

static _NT_factory g_factory = {
    .guid = NT_MULTICHAR('C','s','A','G'),
    .name = "Cessna AWG",
    .description = "CV-driven arbitrary waveform generator",
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
    .midiMessage  = nullptr,
    .tags = kNT_tagUtility,
    .hasCustomUi = hasCustomUi,
    .customUi    = customUi,
    .setupUi     = nullptr,
    .serialise   = nullptr,
    .deserialise = nullptr,
    .midiSysEx   = nullptr,
    .parameterUiPrefix = nullptr,
    .parameterString   = nullptr,
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

