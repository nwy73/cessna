/*
 * Cessna AWG (v2) — CV-driven Arbitrary Waveform Generator
 * ---------------------------------------------------------
 * Wavetable oscillator driven by 1V/oct pitch CV. 0V = C4.
 * 16 manually programmable amplitude bins define waveform shape.
 *
 * Modulation modes (mutually exclusive):
 *   Ripple   — sine wave travelling across bins (per-sample)
 *   Kink     — a peak or notch scanning across bins (timed)
 *   Scramble — randomly swap adjacent bin pairs (timed)
 *   Seq      — pattern rotating through bins (timed)
 *   Drunk    — whole waveform coherently random-walks (timed)
 *
 * Morph: crossfade between two stored waveform snapshots (A/B).
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
static constexpr float kC4Hz             = 261.63f;

// ----------------------------
// Helpers
// ----------------------------
static inline float clampf(float x, float lo, float hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}
static inline float lerpf(float a, float b, float t) {
    return a + (b - a) * t;
}
static inline float lcgRand(uint32_t &state) {
    state = state * 1664525u + 1013904223u;
    return ((float)(state >> 1) / (float)0x7fffffff) - 1.0f; // -1..1
}
static inline float lcgRand01(uint32_t &state) {
    state = state * 1664525u + 1013904223u;
    return (float)(state >> 1) / (float)0x7fffffff; // 0..1
}

// ----------------------------
// Forward declarations
// ----------------------------
static inline float evalLinear(const float *bins, float phase01);
static inline float evalStepped(const float *bins, float phase01);

// ----------------------------
// Cubic Hermite interpolator
// ----------------------------
struct PeriodicHermite16 {
    float m[kBins] = {0};

    void computeSlopes(const float *p) {
        float d[kBins];
        for (int i=0;i<kBins;i++)
            d[i] = p[(i+1)&15] - p[i];
        for (int i=0;i<kBins;i++)
            m[i] = 0.5f * (d[(i-1)&15] + d[i]);
        for (int i=0;i<kBins;i++) {
            float dp = d[i];
            if (dp == 0.0f) { m[i] = 0.0f; m[(i+1)&15] = 0.0f; continue; }
            if ((m[i] / dp) < 0.0f)        m[i]        = 0.0f;
            if ((m[(i+1)&15] / dp) < 0.0f) m[(i+1)&15] = 0.0f;
        }
        for (int i=0;i<kBins;i++) {
            float dp = d[i];
            if (dp == 0.0f) continue;
            float r1 = m[i] / dp, r2 = m[(i+1)&15] / dp;
            float ss = r1*r1 + r2*r2;
            if (ss > 9.0f) { float t = 3.0f/sqrtf(ss); m[i]*=t; m[(i+1)&15]*=t; }
        }
    }

    inline float eval(const float *p, float phase01) const {
        float x = phase01 * (float)kBins;
        int i = (int)floorf(x) & 15;
        float t = x - floorf(x);
        int i1 = (i+1) & 15;
        float t2=t*t, t3=t2*t;
        return (2*t3-3*t2+1)*p[i] + (t3-2*t2+t)*m[i] + (-2*t3+3*t2)*p[i1] + (t3-t2)*m[i1];
    }
};

// ----------------------------
// DTC
// ----------------------------
struct _CessnaAwg_DTC {
    float sampleRate;

    // Waveform bins — direct amplitude values, ±5V
    float B[kBins];
    float S[kBins];

    float depthD;
    float curve;
    bool  stepped;

    float phase;
    float phaseInc;

    PeriodicHermite16 herm;

    // --- Modulation ---
    int   modMode;      // 0=Off 1=Ripple 2=Kink 3=Scramble 4=Seq 5=Drunk

    // Ripple (per-sample LFO)
    float rippleAmt;
    float rippleRate;
    float rippleLfoPhase;
    float rippleLfoPhaseInc;

    // Kink (timed scan)
    float kinkAmt;
    float kinkWidth;    // 0..1 fraction of bins
    bool  kinkIsNotch;  // false=peak, true=notch
    int   kinkPos;      // current bin position (0..kBins-1)
    int   kinkSampleCount;
    int   kinkSamplePeriod;

    // Scramble (timed)
    float scrambleAmt;
    int   scrambleSampleCount;
    int   scrambleSamplePeriod;
    float scrambledBins[kBins]; // working copy that gets scrambled

    // Seq (timed)
    float seqAmt;
    float seqSmooth;
    int   seqPattern;   // 0=Ramp 1=Triangle 2=Alternating
    float seqTarget[kBins];
    float seqCurrent[kBins];
    int   seqSampleCount;
    int   seqSamplePeriod;
    uint32_t seqRng;

    // Drunk (timed coherent walk)
    float drunkAmt;
    float drunkOffset;  // single value added to all bins
    float drunkStep;    // 0..1 max step size per advance
    int   drunkSampleCount;
    int   drunkSamplePeriod;
    uint32_t drunkRng;

    float modX[kBins];

    // Morph slots
    float slotA[kBins];
    float slotB[kBins];
    bool  saveA;
    bool  saveB;
    float morph;

    // Display
    float dispWave[kDisplayPoints];

    // UI
    bool  needsGrayOutUpdate;
    int   algorithmIdx;

    // Compute sample period from rate param (1..32, higher=faster)
    // Maps exponentially: 1 → ~2s, 32 → ~25ms
    static int rateToPeriod(int param, float sr) {
        float pct = (float)param / 32.0f;
        float secs = 2.0f * powf(0.0125f / 2.0f, pct);
        int p = (int)(sr * secs);
        return (p < 1) ? 1 : p;
    }

    _CessnaAwg_DTC() {
        sampleRate = kDefaultSampleRate;
        depthD  = 1.0f;
        curve   = 0.0f;
        stepped = false;
        phase   = 0.0f;
        phaseInc = kC4Hz / kDefaultSampleRate;

        modMode = 0;

        // Ripple
        rippleAmt        = 0.5f;
        rippleRate       = 0.5f;
        rippleLfoPhase   = 0.0f;
        rippleLfoPhaseInc = rippleRate / kDefaultSampleRate;

        // Kink
        kinkAmt          = 0.5f;
        kinkWidth        = 0.125f; // 2 bins wide
        kinkIsNotch      = false;
        kinkPos          = 0;
        kinkSampleCount  = 0;
        kinkSamplePeriod = rateToPeriod(16, kDefaultSampleRate);

        // Scramble
        scrambleAmt          = 0.5f;
        scrambleSampleCount  = 0;
        scrambleSamplePeriod = rateToPeriod(16, kDefaultSampleRate);

        // Seq
        seqAmt          = 0.5f;
        seqSmooth       = 0.5f;
        seqPattern      = 0;
        seqSampleCount  = 0;
        seqSamplePeriod = rateToPeriod(16, kDefaultSampleRate);
        seqRng          = 0x12345678u;

        // Drunk
        drunkAmt          = 0.5f;
        drunkOffset       = 0.0f;
        drunkStep         = 0.05f;
        drunkSampleCount  = 0;
        drunkSamplePeriod = rateToPeriod(8, kDefaultSampleRate);
        drunkRng          = 0x87654321u;

        for (int i=0;i<kBins;i++) {
            B[i] = S[i] = modX[i] = 0.0f;
            slotA[i] = slotB[i] = 0.0f;
            scrambledBins[i] = 0.0f;
        }

        // Seq: initialise ramp pattern, copy to current
        for (int i=0;i<kBins;i++) {
            seqTarget[i]  = ((float)i / (float)(kBins-1)) * 2.0f - 1.0f;
            seqCurrent[i] = seqTarget[i];
        }

        saveA = saveB = false;
        morph = 0.0f;

        for (int i=0;i<kDisplayPoints;i++) dispWave[i] = 0.0f;

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

    // Main page
    kParamDepth,
    kParamCurve,
    kParamStepped,

    // Waveform page (bins)
    kParamB1,
    kParamB16 = kParamB1 + 15,

    // Mod page
    kParamModMode,

    // Ripple page
    kParamRippleAmt,
    kParamRippleRate,

    // Kink page
    kParamKinkAmt,
    kParamKinkRate,
    kParamKinkWidth,
    kParamKinkShape,    // 0=Peak 1=Notch

    // Scramble page
    kParamScrambleAmt,
    kParamScrambleRate,

    // Seq page
    kParamSeqAmt,
    kParamSeqRate,
    kParamSeqSmooth,
    kParamSeqPattern,

    // Drunk page
    kParamDrunkAmt,
    kParamDrunkRate,
    kParamDrunkStep,

    // Morph page
    kParamSaveA,
    kParamSaveB,
    kParamMorph,
};
static constexpr int kNumParams = kParamMorph + 1;

// ----------------------------
// String tables
// ----------------------------
static const char* modModeStrings[]  = {"Off","Ripple","Kink","Scramble","Seq","Drunk",nullptr};
static const char* seqPatStrings[]   = {"Ramp","Triangle","Alternating",nullptr};
static const char* kinkShapeStrings[]= {"Peak","Notch",nullptr};
static const char* onOffStrings[]    = {"Off","On",nullptr};

// ----------------------------
// Parameter descriptors
// ----------------------------
static const _NT_parameter g_parameters[kNumParams] = {
    NT_PARAMETER_CV_INPUT("CV in", 0, 1)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("out", 0, 13)

    // Main
    { .name="Level",    .min=0,  .max=100,.def=100,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Smoothing",.min=0,  .max=100,.def=0,  .unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Stepped",  .min=0,  .max=1,  .def=0,  .unit=kNT_unitEnum,  .scaling=kNT_scalingNone,.enumStrings=onOffStrings },

    // Bins 1–16
    { .name="Bin 1", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 2", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 3", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 4", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 5", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 6", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 7", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 8", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 9", .min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 10",.min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 11",.min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 12",.min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 13",.min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 14",.min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 15",.min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Bin 16",.min=-100,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },

    // Mod selector
    { .name="Mod Mode",.min=0,.max=5,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=modModeStrings },

    // Ripple
    { .name="Amt",  .min=0,.max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Rate", .min=1,.max=1000,.def=50,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },

    // Kink
    { .name="Amt",  .min=0,.max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Rate", .min=1,.max=32,.def=16,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Width",.min=1,.max=8,.def=2,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Shape",.min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=kinkShapeStrings },

    // Scramble
    { .name="Amt",  .min=0,.max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Rate", .min=1,.max=32,.def=16,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },

    // Seq
    { .name="Amt",    .min=0,.max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Rate",   .min=1,.max=32,.def=16,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Smooth", .min=0,.max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Pattern",.min=0,.max=2,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=seqPatStrings },

    // Drunk
    { .name="Amt",  .min=0,.max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Rate", .min=1,.max=32,.def=8,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Step", .min=1,.max=100,.def=10,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },

    // Morph
    { .name="Save A",.min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings },
    { .name="Save B",.min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings },
    { .name="Morph", .min=0,.max=100,.def=0,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
};

// ----------------------------
// Page index arrays
// ----------------------------
static const uint8_t g_pageRoutingIdx[] = {
    (uint8_t)kParamCVIn,
    (uint8_t)kParamOut1,(uint8_t)kParamOut1Mode,
};
static const uint8_t g_pageMainIdx[] = {
    (uint8_t)kParamDepth,(uint8_t)kParamCurve,(uint8_t)kParamStepped,
};
static const uint8_t g_pageBinsIdx[] = {
    (uint8_t)kParamB1,      (uint8_t)(kParamB1+1),  (uint8_t)(kParamB1+2),  (uint8_t)(kParamB1+3),
    (uint8_t)(kParamB1+4),  (uint8_t)(kParamB1+5),  (uint8_t)(kParamB1+6),  (uint8_t)(kParamB1+7),
    (uint8_t)(kParamB1+8),  (uint8_t)(kParamB1+9),  (uint8_t)(kParamB1+10), (uint8_t)(kParamB1+11),
    (uint8_t)(kParamB1+12), (uint8_t)(kParamB1+13), (uint8_t)(kParamB1+14), (uint8_t)(kParamB1+15),
};
static const uint8_t g_pageModIdx[]      = { (uint8_t)kParamModMode };
static const uint8_t g_pageRippleIdx[]   = { (uint8_t)kParamRippleAmt,  (uint8_t)kParamRippleRate };
static const uint8_t g_pageKinkIdx[]     = { (uint8_t)kParamKinkAmt, (uint8_t)kParamKinkRate, (uint8_t)kParamKinkWidth, (uint8_t)kParamKinkShape };
static const uint8_t g_pageScrambleIdx[] = { (uint8_t)kParamScrambleAmt, (uint8_t)kParamScrambleRate };
static const uint8_t g_pageSeqIdx[]      = { (uint8_t)kParamSeqAmt, (uint8_t)kParamSeqRate, (uint8_t)kParamSeqSmooth, (uint8_t)kParamSeqPattern };
static const uint8_t g_pageDrunkIdx[]    = { (uint8_t)kParamDrunkAmt, (uint8_t)kParamDrunkRate, (uint8_t)kParamDrunkStep };
static const uint8_t g_pageMorphIdx[]    = { (uint8_t)kParamSaveA, (uint8_t)kParamSaveB, (uint8_t)kParamMorph };

static const _NT_parameterPage g_pages_arr[] = {
    { .name="Routing", .numParams=(uint8_t)ARRAY_SIZE(g_pageRoutingIdx), .group=0,.unused={0,0},.params=g_pageRoutingIdx },
    { .name="Main",    .numParams=(uint8_t)ARRAY_SIZE(g_pageMainIdx),    .group=0,.unused={0,0},.params=g_pageMainIdx },
    { .name="Waveform",.numParams=(uint8_t)ARRAY_SIZE(g_pageBinsIdx),    .group=1,.unused={0,0},.params=g_pageBinsIdx },
    { .name="Mod",     .numParams=(uint8_t)ARRAY_SIZE(g_pageModIdx),     .group=0,.unused={0,0},.params=g_pageModIdx },
    { .name="Ripple",  .numParams=(uint8_t)ARRAY_SIZE(g_pageRippleIdx),  .group=0,.unused={0,0},.params=g_pageRippleIdx },
    { .name="Kink",    .numParams=(uint8_t)ARRAY_SIZE(g_pageKinkIdx),    .group=0,.unused={0,0},.params=g_pageKinkIdx },
    { .name="Scramble",.numParams=(uint8_t)ARRAY_SIZE(g_pageScrambleIdx),.group=0,.unused={0,0},.params=g_pageScrambleIdx },
    { .name="Seq",     .numParams=(uint8_t)ARRAY_SIZE(g_pageSeqIdx),     .group=0,.unused={0,0},.params=g_pageSeqIdx },
    { .name="Drunk",   .numParams=(uint8_t)ARRAY_SIZE(g_pageDrunkIdx),   .group=0,.unused={0,0},.params=g_pageDrunkIdx },
    { .name="Morph",   .numParams=(uint8_t)ARRAY_SIZE(g_pageMorphIdx),   .group=0,.unused={0,0},.params=g_pageMorphIdx },
};

static const _NT_parameterPages g_pages = {
    .numPages = ARRAY_SIZE(g_pages_arr),
    .pages    = g_pages_arr,
};

// ----------------------------
// Helpers
// ----------------------------
#define getPct01(p)     ((float)self->v[(p)] * 0.01f)
#define getPctSigned(p) ((float)self->v[(p)] * 0.01f * 5.0f)

static inline float evalLinear(const float *bins, float phase01) {
    float x = phase01 * (float)kBins;
    int   i = (int)floorf(x) & 15;
    float t = x - floorf(x);
    return lerpf(bins[i], bins[(i+1)&15], t);
}
static inline float evalStepped(const float *bins, float phase01) {
    return bins[(int)floorf(phase01 * (float)kBins) & 15];
}

// ----------------------------
// combineS
// ----------------------------
static void combineS(struct _CessnaAwg_DTC *d) {
    if (!d->saveA) for (int i=0;i<kBins;i++) d->slotA[i] = d->B[i];
    if (!d->saveB) for (int i=0;i<kBins;i++) d->slotB[i] = d->B[i];
    float D = d->depthD;
    for (int i=0;i<kBins;i++) {
        float base = lerpf(d->slotA[i], d->slotB[i], d->morph);
        d->S[i] = D * tanhf((base + d->modX[i] * 5.0f) / 5.0f) * 5.0f;
    }
}

// ----------------------------
// Modulation functions
// ----------------------------

// RIPPLE: sine LFO phase-offset per bin, advances every sample
static void computeRipple(struct _CessnaAwg_DTC *d) {
    d->rippleLfoPhase += d->rippleLfoPhaseInc;
    if (d->rippleLfoPhase >= 1.0f) d->rippleLfoPhase -= 1.0f;
    for (int i=0;i<kBins;i++) {
        float bp = d->rippleLfoPhase + (float)i / (float)kBins;
        if (bp >= 1.0f) bp -= 1.0f;
        d->modX[i] = d->rippleAmt * sinf(2.0f * M_PI_F * bp);
    }
}

// KINK: a peak or notch of configurable width scans across bins
// modX is built fresh each sample from current kinkPos
static void buildKinkModX(struct _CessnaAwg_DTC *d) {
    float polarity = d->kinkIsNotch ? -1.0f : 1.0f;
    int   halfW    = (int)roundf(d->kinkWidth * (float)kBins * 0.5f);
    if (halfW < 1) halfW = 1;
    for (int i=0;i<kBins;i++) d->modX[i] = 0.0f;
    for (int w=-halfW; w<=halfW; w++) {
        int idx = (d->kinkPos + w + kBins * 4) % kBins;
        // Triangular shape: full amplitude at centre, zero at edges
        float env = 1.0f - fabsf((float)w / (float)(halfW + 1));
        d->modX[idx] = polarity * d->kinkAmt * env;
    }
}

static void advanceKink(struct _CessnaAwg_DTC *d) {
    d->kinkPos = (d->kinkPos + 1) % kBins;
}

// SCRAMBLE: on each step, randomly swap some adjacent bin pairs
// Works on a working copy (scrambledBins) that resets to B[] each step
static void advanceScramble(struct _CessnaAwg_DTC *d) {
    // Reset working copy to current base bins
    for (int i=0;i<kBins;i++) d->scrambledBins[i] = d->B[i];
    // Randomly swap adjacent pairs — number of swaps scales with amt
    int swaps = 1 + (int)(d->scrambleAmt * (float)(kBins/2 - 1));
    for (int s=0;s<swaps;s++) {
        uint32_t r = 0;
        lcgRand01(d->seqRng); // reuse seqRng for scramble randomness
        r = d->seqRng;
        int idx = (int)(lcgRand01(r) * (float)(kBins-1));
        if (idx >= kBins-1) idx = kBins-2;
        float tmp = d->scrambledBins[idx];
        d->scrambledBins[idx]   = d->scrambledBins[idx+1];
        d->scrambledBins[idx+1] = tmp;
    }
}

static void buildScrambleModX(struct _CessnaAwg_DTC *d) {
    // modX = difference between scrambled and original, scaled by amt
    for (int i=0;i<kBins;i++)
        d->modX[i] = (d->scrambledBins[i] - d->B[i]) / 5.0f; // normalise to ±1
}

// SEQ: a pattern rotates through the bins
static void buildSeqPattern(struct _CessnaAwg_DTC *d) {
    switch (d->seqPattern) {
        case 0: // Ramp
            for (int i=0;i<kBins;i++)
                d->seqTarget[i] = ((float)i / (float)(kBins-1)) * 2.0f - 1.0f;
            break;
        case 1: // Triangle
            for (int i=0;i<kBins;i++) {
                float t = (float)i / (float)kBins;
                d->seqTarget[i] = (t < 0.5f) ? (t*4.0f-1.0f) : (3.0f-t*4.0f);
            }
            break;
        case 2: // Alternating
            for (int i=0;i<kBins;i++)
                d->seqTarget[i] = (i & 1) ? -1.0f : 1.0f;
            break;
    }
    // Copy to current immediately to avoid initial slew from zero
    for (int i=0;i<kBins;i++) d->seqCurrent[i] = d->seqTarget[i];
}

static void advanceSeq(struct _CessnaAwg_DTC *d) {
    // Rotate target array by one position
    float tmp = d->seqTarget[kBins-1];
    for (int i=kBins-1;i>0;i--) d->seqTarget[i] = d->seqTarget[i-1];
    d->seqTarget[0] = tmp;
}

static void computeSeqModX(struct _CessnaAwg_DTC *d) {
    // Slew seqCurrent toward seqTarget each sample
    float slew = 1.0f - d->seqSmooth * 0.999f;
    for (int i=0;i<kBins;i++) {
        d->seqCurrent[i] += slew * (d->seqTarget[i] - d->seqCurrent[i]);
        d->modX[i] = d->seqAmt * d->seqCurrent[i];
    }
}

// DRUNK: single random offset applied coherently to all bins
// Preserves waveform shape, animates overall amplitude
static void advanceDrunk(struct _CessnaAwg_DTC *d) {
    float step = lcgRand(d->drunkRng) * d->drunkStep;
    d->drunkOffset = clampf(d->drunkOffset + step, -1.0f, 1.0f);
}

static void buildDrunkModX(struct _CessnaAwg_DTC *d) {
    float v = d->drunkAmt * d->drunkOffset;
    for (int i=0;i<kBins;i++) d->modX[i] = v;
}

// Master computeModX dispatcher
static void computeModX(struct _CessnaAwg_DTC *d) {
    switch (d->modMode) {
        case 0: for (int i=0;i<kBins;i++) d->modX[i] = 0.0f; break;
        case 1: computeRipple(d);     break;
        case 2: buildKinkModX(d);     break;
        case 3: buildScrambleModX(d); break;
        case 4: computeSeqModX(d);    break;
        case 5: buildDrunkModX(d);    break;
    }
}

// ----------------------------
// parameterChanged()
// ----------------------------
static void parameterChanged(_NT_algorithm *base, int p) {
    auto *self = (_CessnaAwg*)base;
    auto *d = self->dtc;
    if (!d) return;

    switch (p) {
        case kParamDepth:   d->depthD  = getPct01(kParamDepth); break;
        case kParamCurve:   d->curve   = getPct01(kParamCurve); break;
        case kParamStepped: d->stepped = (self->v[kParamStepped] != 0); break;

        case kParamModMode: d->modMode = self->v[kParamModMode]; break;

        // Ripple
        case kParamRippleAmt:  d->rippleAmt = getPct01(kParamRippleAmt); break;
        case kParamRippleRate: {
            float hz = (float)self->v[kParamRippleRate] * 0.01f;
            d->rippleRate = hz;
            d->rippleLfoPhaseInc = hz / d->sampleRate;
            break;
        }

        // Kink
        case kParamKinkAmt:   d->kinkAmt     = getPct01(kParamKinkAmt); break;
        case kParamKinkRate:  d->kinkSamplePeriod  = _CessnaAwg_DTC::rateToPeriod(self->v[kParamKinkRate], d->sampleRate); break;
        case kParamKinkWidth: d->kinkWidth   = (float)self->v[kParamKinkWidth] / (float)kBins; break;
        case kParamKinkShape: d->kinkIsNotch = (self->v[kParamKinkShape] != 0); break;

        // Scramble
        case kParamScrambleAmt:  d->scrambleAmt          = getPct01(kParamScrambleAmt); break;
        case kParamScrambleRate: d->scrambleSamplePeriod = _CessnaAwg_DTC::rateToPeriod(self->v[kParamScrambleRate], d->sampleRate); break;

        // Seq
        case kParamSeqAmt:     d->seqAmt     = getPct01(kParamSeqAmt); break;
        case kParamSeqRate:    d->seqSamplePeriod = _CessnaAwg_DTC::rateToPeriod(self->v[kParamSeqRate], d->sampleRate); break;
        case kParamSeqSmooth:  d->seqSmooth  = getPct01(kParamSeqSmooth); break;
        case kParamSeqPattern: {
            d->seqPattern = self->v[kParamSeqPattern];
            buildSeqPattern(d);
            break;
        }

        // Drunk
        case kParamDrunkAmt:  d->drunkAmt          = getPct01(kParamDrunkAmt); break;
        case kParamDrunkRate: d->drunkSamplePeriod = _CessnaAwg_DTC::rateToPeriod(self->v[kParamDrunkRate], d->sampleRate); break;
        case kParamDrunkStep: d->drunkStep          = getPct01(kParamDrunkStep); break;

        // Morph
        case kParamSaveA: {
            d->saveA = (self->v[kParamSaveA] != 0);
            if (d->saveA) for (int i=0;i<kBins;i++) d->slotA[i] = d->B[i];
            break;
        }
        case kParamSaveB: {
            d->saveB = (self->v[kParamSaveB] != 0);
            if (d->saveB) for (int i=0;i<kBins;i++) d->slotB[i] = d->B[i];
            break;
        }
        case kParamMorph: d->morph = getPct01(kParamMorph); break;

        default:
            if (p >= kParamB1 && p <= kParamB16) {
                d->B[p - kParamB1] = getPctSigned(p);
                // Reset scramble working copy when bins change
                for (int i=0;i<kBins;i++) d->scrambledBins[i] = d->B[i];
            }
            break;
    }
}

// ----------------------------
// updateDisplay
// ----------------------------
static void updateDisplay(struct _CessnaAwg_DTC *d) {
    for (int k=0;k<kDisplayPoints;k++) {
        float ph = (float)k / (float)(kDisplayPoints-1);
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

    const int cvBus  = self->v[kParamCVIn];
    const int outBus = self->v[kParamOut1];
    const int outMode= self->v[kParamOut1Mode];

    auto busPtr = [&](int busVal) -> float* {
        if (!busFrames || busVal <= 0) return nullptr;
        int idx = busVal - 1;
        if (idx < 0 || idx >= 28) return nullptr;
        return busFrames + idx * numFrames;
    };

    float *cv  = busPtr(cvBus);
    float *out = busPtr(outBus);

    for (int n=0;n<numFrames;n++) {
        // Pitch
        float cvVal = cv ? cv[n] : 0.0f;
        d->phaseInc = kC4Hz * powf(2.0f, cvVal) / d->sampleRate;

        // Timed mod advances (Kink, Scramble, Seq, Drunk)
        if (d->modMode == 2) { // Kink
            d->kinkSampleCount++;
            if (d->kinkSampleCount >= d->kinkSamplePeriod) {
                d->kinkSampleCount = 0;
                advanceKink(d);
            }
        } else if (d->modMode == 3) { // Scramble
            d->scrambleSampleCount++;
            if (d->scrambleSampleCount >= d->scrambleSamplePeriod) {
                d->scrambleSampleCount = 0;
                advanceScramble(d);
            }
        } else if (d->modMode == 4) { // Seq
            d->seqSampleCount++;
            if (d->seqSampleCount >= d->seqSamplePeriod) {
                d->seqSampleCount = 0;
                advanceSeq(d);
            }
        } else if (d->modMode == 5) { // Drunk
            d->drunkSampleCount++;
            if (d->drunkSampleCount >= d->drunkSamplePeriod) {
                d->drunkSampleCount = 0;
                advanceDrunk(d);
            }
        }

        // Compute modX for this sample
        computeModX(d);

        // Combine bins + mod + morph → S[]
        combineS(d);
        if (!d->stepped && d->curve > 0.0001f)
            d->herm.computeSlopes(d->S);

        // Evaluate wavetable
        float y;
        if (d->stepped) {
            y = evalStepped(d->S, d->phase);
        } else {
            float lin = evalLinear(d->S, d->phase);
            float cur = (d->curve > 0.0001f) ? d->herm.eval(d->S, d->phase) : lin;
            y = lerpf(lin, cur, d->curve);
        }

        if (out) {
            if (outMode == 0) out[n] += y;
            else              out[n]  = y;
        }

        d->phase += d->phaseInc;
        if (d->phase >= 1.0f) d->phase -= 1.0f;
    }

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
    const int top = 11, bottom = H-2;
    const int plotH = bottom - top;

    auto mapY = [&](float v)->int {
        float vv = clampf(v, -5.0f, 5.0f);
        int y = top + (int)roundf(0.5f*(1.0f-(vv/5.0f)) * (float)plotH);
        return clampf((float)y, (float)top, (float)(H-1));
    };

    for (int n=1;n<kDisplayPoints;n++) {
        int x0 = (int)roundf(((float)(n-1)/(float)(kDisplayPoints-1))*(float)(W-1));
        int x1 = (int)roundf(((float)n    /(float)(kDisplayPoints-1))*(float)(W-1));
        int y0 = mapY(d->dispWave[n-1]);
        int y1 = mapY(d->dispWave[n]);
        NT_drawShapeI(kNT_line, x0, y0, x1, y1, 15);
        NT_drawShapeI(kNT_line, x0, y0, x1, y1, 15);
    }

    // Bin tick marks
    for (int i=0;i<kBins;i++) {
        float ph = ((float)i+0.5f)/(float)kBins;
        int x = (int)roundf(ph*(float)(W-1));
        NT_drawShapeI(kNT_line, x, H-4, x, H-1, 6);
    }

    return false;
}

// ----------------------------
// NT registration
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
    self->parameters     = g_parameters;
    self->parameterPages = &g_pages;

    float sr = (float)NT_globals.sampleRate;
    if (sr < 8000.0f || sr > 192000.0f) sr = kDefaultSampleRate;
    d->sampleRate         = sr;
    d->rippleLfoPhaseInc  = d->rippleRate / sr;
    d->kinkSamplePeriod   = _CessnaAwg_DTC::rateToPeriod(16, sr);
    d->scrambleSamplePeriod = _CessnaAwg_DTC::rateToPeriod(16, sr);
    d->seqSamplePeriod    = _CessnaAwg_DTC::rateToPeriod(16, sr);
    d->drunkSamplePeriod  = _CessnaAwg_DTC::rateToPeriod(8, sr);
    d->phaseInc           = kC4Hz / sr;

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
    .initialise   = initialise,
    .calculateRequirements = calculateRequirements,
    .construct    = construct,
    .parameterChanged = parameterChanged,
    .step         = step,
    .draw         = draw,
    .midiRealtime = nullptr,
    .midiMessage  = nullptr,
    .tags         = kNT_tagUtility,
    .hasCustomUi  = hasCustomUi,
    .customUi     = customUi,
    .setupUi      = nullptr,
    .serialise    = nullptr,
    .deserialise  = nullptr,
    .midiSysEx    = nullptr,
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
