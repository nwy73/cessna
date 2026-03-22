/*
 * Cessna AWG (v3) — CV-driven Arbitrary Waveform Generator
 * ---------------------------------------------------------
 * Wavetable oscillator driven by 1V/oct pitch CV. 0V = C4.
 * 16 manually programmable amplitude bins define waveform shape.
 *
 * Modulation modes:
 *   Wave     — a shaped perturbation travels across the bins.
 *              Shape:     Sine / Triangle / Pulse
 *              Direction: Forward / Backward / Alternating / Random / Drunk
 *              Width:     how wide the perturbation is (bins)
 *              Range:     how far it travels (bins)
 *              Rate:      speed
 *              Amt:       depth
 *   Scramble — randomly swap adjacent bin pairs each step
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
    // Returns -1..1
    state = state * 1664525u + 1013904223u;
    return ((float)(state >> 1) / (float)0x7fffffff) - 1.0f;
}
static inline float lcgRand01(uint32_t &state) {
    // Returns 0..1
    state = state * 1664525u + 1013904223u;
    return (float)(state >> 1) / (float)0x7fffffff;
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
        return (2*t3-3*t2+1)*p[i] + (t3-2*t2+t)*m[i] +
               (-2*t3+3*t2)*p[i1] + (t3-t2)*m[i1];
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

    // --- Modulation mode ---
    // 0=Off  1=Ripple  2=Wave  3=Scramble
    int   modMode;

    // --- Ripple mod ---
    float rippleAmt;
    float rippleLfoPhase;
    float rippleLfoPhaseInc;

    // --- Wave mod ---
    float waveAmt;       // 0..1 depth
    int   waveShape;     // 0=Sine 1=Triangle 2=Pulse
    int   waveWidth;     // 1..16 bins, width of perturbation
    int   waveRange;     // 1..16 bins, travel distance
    int   waveDirection; // 0=Fwd 1=Back 2=Alt 3=Random 4=Drunk
    int   waveSampleCount;
    int   waveSamplePeriod;

    // Wave position state
    float wavePos;       // current centre position (float, 0..kBins-1)
    float waveDrunkPos;  // for Drunk direction: current pos as float
    float waveDrunkStep; // max step per advance for Drunk
    int   waveAltDir;    // +1 or -1 for Alternating direction
    int   waveRangeStart;// leftmost bin of travel range (set on first advance)
    bool  waveRangeInit; // false until first advance initialises range

    float modX[kBins];

    // --- Scramble mod ---
    float scrambleAmt;
    int   scrambleSampleCount;
    int   scrambleSamplePeriod;
    uint32_t scrambleRng;

    // --- Morph slots ---
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

    // Map rate param (1..32, higher=faster) to sample period
    // 1 → ~2s, 32 → ~25ms (exponential)
    static int rateToPeriod(int param, float sr) {
        float pct  = (float)param / 32.0f;
        float secs = 2.0f * powf(0.0125f / 2.0f, pct);
        int   p    = (int)(sr * secs);
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

        // Ripple defaults
        rippleAmt          = 0.5f;
        rippleLfoPhase     = 0.0f;
        rippleLfoPhaseInc  = 0.5f / kDefaultSampleRate; // 0.5Hz default

        // Wave defaults
        waveAmt         = 0.5f;
        waveShape       = 0;    // Sine
        waveWidth       = 2;    // 2 bins wide
        waveRange       = kBins; // full range
        waveDirection   = 0;    // Forward
        waveSampleCount = 0;
        waveSamplePeriod = rateToPeriod(16, kDefaultSampleRate);
        wavePos         = 0.0f;
        waveDrunkPos    = 0.0f;
        waveDrunkStep   = 1.0f;
        waveAltDir      = 1;
        waveRangeStart  = 0;
        waveRangeInit   = false;

        // Scramble defaults
        scrambleAmt          = 0.5f;
        scrambleSampleCount  = 0;
        scrambleSamplePeriod = rateToPeriod(16, kDefaultSampleRate);
        scrambleRng          = 0xABCD1234u;

        for (int i=0;i<kBins;i++) {
            B[i] = S[i] = modX[i] = 0.0f;
            slotA[i] = slotB[i] = 0.0f;
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

    // Waveform page (bins 1-16)
    kParamB1,
    kParamB16 = kParamB1 + 15,

    // Mod selector page
    kParamModMode,

    // Ripple page
    kParamRippleAmt,
    kParamRippleRate,

    // Wave page
    kParamWaveAmt,
    kParamWaveRate,
    kParamWaveShape,
    kParamWaveWidth,
    kParamWaveRange,
    kParamWaveDirection,

    // Scramble page
    kParamScrambleAmt,
    kParamScrambleRate,

    // Morph page
    kParamSaveA,
    kParamSaveB,
    kParamMorph,
};
static constexpr int kNumParams = kParamMorph + 1;

// ----------------------------
// String tables
// ----------------------------
static const char* modModeStrings[]  = {"Off","Ripple","Wave","Scramble",nullptr};
static const char* waveShapeStrings[]= {"Triangle","Pulse",nullptr};
static const char* waveDirStrings[]  = {"Forward","Backward","Alternating","Random","Drunk",nullptr};
static const char* onOffStrings[]    = {"Off","On",nullptr};

// ----------------------------
// Parameter descriptors
// (must be in enum index order)
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
    { .name="Mod Mode",.min=0,.max=3,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=modModeStrings },

    // Ripple
    { .name="Amt",  .min=0,.max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Rate", .min=1,.max=1000,.def=50,.unit=kNT_unitNone,  .scaling=kNT_scalingNone,.enumStrings=nullptr },

    // Wave
    { .name="Amt",      .min=0,.max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Rate",     .min=1,.max=32, .def=16,.unit=kNT_unitNone,  .scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Shape",    .min=0,.max=1,  .def=0, .unit=kNT_unitEnum,  .scaling=kNT_scalingNone,.enumStrings=waveShapeStrings },
    { .name="Width",    .min=1,.max=16, .def=2, .unit=kNT_unitNone,  .scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Range",    .min=1,.max=16, .def=16,.unit=kNT_unitNone,  .scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Direction",.min=0,.max=4,  .def=0, .unit=kNT_unitEnum,  .scaling=kNT_scalingNone,.enumStrings=waveDirStrings },

    // Scramble
    { .name="Amt",  .min=0,.max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Rate", .min=1,.max=32, .def=16,.unit=kNT_unitNone,  .scaling=kNT_scalingNone,.enumStrings=nullptr },

    // Morph
    { .name="Save A",.min=0,.max=1,  .def=0,.unit=kNT_unitEnum,   .scaling=kNT_scalingNone,.enumStrings=onOffStrings },
    { .name="Save B",.min=0,.max=1,  .def=0,.unit=kNT_unitEnum,   .scaling=kNT_scalingNone,.enumStrings=onOffStrings },
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
    (uint8_t)(kParamB1+0),(uint8_t)(kParamB1+1),(uint8_t)(kParamB1+2),(uint8_t)(kParamB1+3),
    (uint8_t)(kParamB1+4),(uint8_t)(kParamB1+5),(uint8_t)(kParamB1+6),(uint8_t)(kParamB1+7),
    (uint8_t)(kParamB1+8),(uint8_t)(kParamB1+9),(uint8_t)(kParamB1+10),(uint8_t)(kParamB1+11),
    (uint8_t)(kParamB1+12),(uint8_t)(kParamB1+13),(uint8_t)(kParamB1+14),(uint8_t)(kParamB1+15),
};
static const uint8_t g_pageModIdx[]      = { (uint8_t)kParamModMode };
static const uint8_t g_pageRippleIdx[]   = { (uint8_t)kParamRippleAmt, (uint8_t)kParamRippleRate };
static const uint8_t g_pageWaveIdx[]     = {
    (uint8_t)kParamWaveAmt,(uint8_t)kParamWaveRate,(uint8_t)kParamWaveShape,
    (uint8_t)kParamWaveWidth,(uint8_t)kParamWaveRange,(uint8_t)kParamWaveDirection,
};
static const uint8_t g_pageScrambleIdx[] = { (uint8_t)kParamScrambleAmt,(uint8_t)kParamScrambleRate };
static const uint8_t g_pageMorphIdx[]    = { (uint8_t)kParamSaveA,(uint8_t)kParamSaveB,(uint8_t)kParamMorph };

static const _NT_parameterPage g_pages_arr[] = {
    { .name="Routing", .numParams=(uint8_t)ARRAY_SIZE(g_pageRoutingIdx),.group=0,.unused={0,0},.params=g_pageRoutingIdx },
    { .name="Main",    .numParams=(uint8_t)ARRAY_SIZE(g_pageMainIdx),   .group=0,.unused={0,0},.params=g_pageMainIdx },
    { .name="Waveform",.numParams=(uint8_t)ARRAY_SIZE(g_pageBinsIdx),   .group=1,.unused={0,0},.params=g_pageBinsIdx },
    { .name="Mod",     .numParams=(uint8_t)ARRAY_SIZE(g_pageModIdx),    .group=0,.unused={0,0},.params=g_pageModIdx },
    { .name="Ripple",  .numParams=(uint8_t)ARRAY_SIZE(g_pageRippleIdx), .group=0,.unused={0,0},.params=g_pageRippleIdx },
    { .name="Wave",    .numParams=(uint8_t)ARRAY_SIZE(g_pageWaveIdx),   .group=0,.unused={0,0},.params=g_pageWaveIdx },
    { .name="Scramble",.numParams=(uint8_t)ARRAY_SIZE(g_pageScrambleIdx),.group=0,.unused={0,0},.params=g_pageScrambleIdx },
    { .name="Morph",   .numParams=(uint8_t)ARRAY_SIZE(g_pageMorphIdx),  .group=0,.unused={0,0},.params=g_pageMorphIdx },
};

static const _NT_parameterPages g_pages = {
    .numPages = ARRAY_SIZE(g_pages_arr),
    .pages    = g_pages_arr,
};

// ----------------------------
// Macros
// ----------------------------
#define getPct01(p)     ((float)self->v[(p)] * 0.01f)
#define getPctSigned(p) ((float)self->v[(p)] * 0.01f * 5.0f)

// ----------------------------
// Wavetable evaluation
// ----------------------------
static inline float evalLinear(const float *bins, float phase01) {
    float x = phase01 * (float)kBins;
    int   i = (int)floorf(x) & 15;
    return lerpf(bins[i], bins[(i+1)&15], x - floorf(x));
}
static inline float evalStepped(const float *bins, float phase01) {
    return bins[(int)floorf(phase01 * (float)kBins) & 15];
}

// ----------------------------
// combineS: slots + morph + modX → S[]
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
// Wave modulation
// ----------------------------

// Build the perturbation shape envelope for a given position
// Returns values in 0..1 (unsigned) — polarity applied by caller
static float waveEnvelope(int shape, int binIdx, int centrePos, int halfWidth) {
    int dist = binIdx - centrePos;
    // Wrap distance to nearest (handles wraparound)
    if (dist >  kBins/2) dist -= kBins;
    if (dist < -kBins/2) dist += kBins;
    int absDist = dist < 0 ? -dist : dist;

    if (absDist > halfWidth) return 0.0f;

    float t = (halfWidth > 0) ? (float)absDist / (float)halfWidth : 0.0f;

    switch (shape) {
        case 0: // Triangle — linear taper
            return 1.0f - t;
        case 1: // Pulse — rectangular, full amplitude across full width
            return 1.0f;
        default:
            return 0.0f;
    }
}

// Advance wave position according to direction
static void advanceWavePosition(struct _CessnaAwg_DTC *d) {
    // Initialise range anchor on first advance
    if (!d->waveRangeInit) {
        d->waveRangeStart = 0;
        d->wavePos        = 0.0f;
        d->waveDrunkPos   = 0.0f;
        d->waveRangeInit  = true;
    }

    float range = (float)(d->waveRange - 1);
    if (range < 0.0f) range = 0.0f;

    switch (d->waveDirection) {
        case 0: // Forward — step +1, wrap within range
            d->wavePos = d->wavePos + 1.0f;
            if (d->wavePos > range) d->wavePos = 0.0f;
            break;
        case 1: // Backward — step -1, wrap within range
            d->wavePos = d->wavePos - 1.0f;
            if (d->wavePos < 0.0f) d->wavePos = range;
            break;
        case 2: // Alternating — pingpong within range
            d->wavePos += (float)d->waveAltDir;
            if (d->wavePos >= range) { d->wavePos = range; d->waveAltDir = -1; }
            if (d->wavePos <= 0.0f) { d->wavePos = 0.0f;  d->waveAltDir =  1; }
            break;
        case 3: { // Random — jump to random position within range
            float r = lcgRand01(d->scrambleRng); // reuse scramble rng
            d->wavePos = r * range;
            break;
        }
        case 4: // Drunk — random walk within range
            d->waveDrunkPos += lcgRand(d->scrambleRng) * d->waveDrunkStep;
            d->waveDrunkPos  = clampf(d->waveDrunkPos, 0.0f, range);
            d->wavePos       = d->waveDrunkPos;
            break;
    }
}

// Build modX from current wave position
static void buildWaveModX(struct _CessnaAwg_DTC *d) {
    int centrePos = ((int)roundf(d->wavePos) + d->waveRangeStart) % kBins;
    int halfWidth = d->waveWidth / 2;
    if (halfWidth < 1) halfWidth = 1;

    for (int i=0;i<kBins;i++) {
        float env = waveEnvelope(d->waveShape, i, centrePos, halfWidth);
        d->modX[i] = d->waveAmt * env;
    }
}

// ----------------------------
// Scramble modulation
// ----------------------------
static void advanceScramble(struct _CessnaAwg_DTC *d) {
    // Number of swaps scales with amt (1 to kBins/2)
    int swaps = 1 + (int)(d->scrambleAmt * (float)(kBins/2 - 1));
    float tmp[kBins];
    for (int i=0;i<kBins;i++) tmp[i] = d->B[i];

    for (int s=0;s<swaps;s++) {
        int idx = (int)(lcgRand01(d->scrambleRng) * (float)(kBins-1));
        if (idx >= kBins-1) idx = kBins-2;
        float t     = tmp[idx];
        tmp[idx]    = tmp[idx+1];
        tmp[idx+1]  = t;
    }
    // modX = difference between scrambled and original, normalised to ±1
    for (int i=0;i<kBins;i++)
        d->modX[i] = (tmp[i] - d->B[i]) / 5.0f;
}

// ----------------------------
// Master computeModX
// ----------------------------
static void computeModX(struct _CessnaAwg_DTC *d) {
    switch (d->modMode) {
        case 0: // Off
            for (int i=0;i<kBins;i++) d->modX[i] = 0.0f;
            break;
        case 1: { // Ripple — continuous sine LFO, all bins phase-offset
            d->rippleLfoPhase += d->rippleLfoPhaseInc;
            if (d->rippleLfoPhase >= 1.0f) d->rippleLfoPhase -= 1.0f;
            for (int i=0;i<kBins;i++) {
                float bp = d->rippleLfoPhase + (float)i / (float)kBins;
                if (bp >= 1.0f) bp -= 1.0f;
                d->modX[i] = d->rippleAmt * sinf(2.0f * M_PI_F * bp);
            }
            break;
        }
        case 2: // Wave — modX already built by advanceWavePosition/buildWaveModX
            buildWaveModX(d);
            break;
        case 3: // Scramble — modX already built by advanceScramble
            break;
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
        // Main
        case kParamDepth:   d->depthD  = getPct01(kParamDepth); break;
        case kParamCurve:   d->curve   = getPct01(kParamCurve); break;
        case kParamStepped: d->stepped = (self->v[kParamStepped] != 0); break;

        // Mod selector
        case kParamModMode:
            d->modMode = self->v[kParamModMode];
            if (d->modMode == 2) { // Wave
                d->wavePos       = 0.0f;
                d->waveDrunkPos  = 0.0f;
                d->waveRangeInit = false;
                d->waveAltDir    = 1;
            }
            break;

        // Ripple
        case kParamRippleAmt:
            d->rippleAmt = getPct01(kParamRippleAmt);
            break;
        case kParamRippleRate: {
            float hz = (float)self->v[kParamRippleRate] * 0.01f;
            d->rippleLfoPhaseInc = hz / d->sampleRate;
            break;
        }

        // Wave
        case kParamWaveAmt:       d->waveAmt       = getPct01(kParamWaveAmt); break;
        case kParamWaveRate:
            d->waveSamplePeriod = _CessnaAwg_DTC::rateToPeriod(self->v[kParamWaveRate], d->sampleRate);
            break;
        case kParamWaveShape:     d->waveShape     = self->v[kParamWaveShape]; break;
        case kParamWaveWidth:     d->waveWidth     = self->v[kParamWaveWidth]; break;
        case kParamWaveRange:
            d->waveRange     = self->v[kParamWaveRange];
            // Clamp current position to new range
            if (d->wavePos >= (float)d->waveRange) d->wavePos = 0.0f;
            if (d->waveDrunkPos >= (float)d->waveRange) d->waveDrunkPos = 0.0f;
            break;
        case kParamWaveDirection:
            d->waveDirection = self->v[kParamWaveDirection];
            // Reset drunk step size — use 1 bin per advance as default
            d->waveDrunkStep = 1.0f;
            break;

        // Scramble
        case kParamScrambleAmt:  d->scrambleAmt          = getPct01(kParamScrambleAmt); break;
        case kParamScrambleRate: d->scrambleSamplePeriod = _CessnaAwg_DTC::rateToPeriod(self->v[kParamScrambleRate], d->sampleRate); break;

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
            if (p >= kParamB1 && p <= kParamB16)
                d->B[p - kParamB1] = getPctSigned(p);
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

        // Timed advances for Wave and Scramble
        if (d->modMode == 2) {
            // Wave: timed step advance
            d->waveSampleCount++;
            if (d->waveSampleCount >= d->waveSamplePeriod) {
                d->waveSampleCount = 0;
                advanceWavePosition(d);
            }
        } else if (d->modMode == 3) {
            // Scramble: timed advance
            d->scrambleSampleCount++;
            if (d->scrambleSampleCount >= d->scrambleSamplePeriod) {
                d->scrambleSampleCount = 0;
                advanceScramble(d);
            }
        }

        // Compute modX
        computeModX(d);

        // Combine → S[]
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

    const int W = 256, H = 64;
    const int top = 11, bottom = H-2, plotH = bottom-top;

    auto mapY = [&](float v)->int {
        float vv = clampf(v, -5.0f, 5.0f);
        int y = top + (int)roundf(0.5f*(1.0f-(vv/5.0f))*(float)plotH);
        if (y < top)  y = top;
        if (y >= H)   y = H-1;
        return y;
    };

    for (int n=1;n<kDisplayPoints;n++) {
        int x0 = (int)roundf(((float)(n-1)/(float)(kDisplayPoints-1))*(float)(W-1));
        int x1 = (int)roundf(((float)n    /(float)(kDisplayPoints-1))*(float)(W-1));
        int y0 = mapY(d->dispWave[n-1]);
        int y1 = mapY(d->dispWave[n]);
        NT_drawShapeI(kNT_line, x0, y0, x1, y1, 15);
        NT_drawShapeI(kNT_line, x0, y0, x1, y1, 15); // double pass = bold
    }

    // Bin tick marks at bottom
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
    d->sampleRate            = sr;
    d->phaseInc              = kC4Hz / sr;
    d->rippleLfoPhaseInc     = 0.5f / sr;
    d->waveSamplePeriod      = _CessnaAwg_DTC::rateToPeriod(16, sr);
    d->scrambleSamplePeriod  = _CessnaAwg_DTC::rateToPeriod(16, sr);

    return self;
}

static uint32_t hasCustomUi(_NT_algorithm *) { return 0; }
static void customUi(_NT_algorithm *, const _NT_uiData &) {}
static void calculateStaticRequirements(_NT_staticRequirements& req) { req.dram = 0; }
static void initialise(_NT_staticMemoryPtrs& ptrs, const _NT_staticRequirements& req) {
    (void)ptrs; (void)req;
}

static _NT_factory g_factory = {
    .guid        = NT_MULTICHAR('C','s','A','G'),
    .name        = "Cessna AWG",
    .description = "CV-driven arbitrary waveform generator",
    .numSpecifications = 0,
    .specifications    = nullptr,
    .calculateStaticRequirements = calculateStaticRequirements,
    .initialise        = initialise,
    .calculateRequirements = calculateRequirements,
    .construct         = construct,
    .parameterChanged  = parameterChanged,
    .step              = step,
    .draw              = draw,
    .midiRealtime      = nullptr,
    .midiMessage       = nullptr,
    .tags              = kNT_tagUtility,
    .hasCustomUi       = hasCustomUi,
    .customUi          = customUi,
    .setupUi           = nullptr,
    .serialise         = nullptr,
    .deserialise       = nullptr,
    .midiSysEx         = nullptr,
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
