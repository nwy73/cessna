/*
 * Cessna AWG (v3) — CV-driven Arbitrary Waveform Generator
 * ---------------------------------------------------------
 * Wavetable oscillator driven by 1V/oct pitch CV. 0V = C4.
 * 16 manually programmable amplitude bins define waveform shape.
 *
 * Modulation modes:
 *   Wave     — a shaped perturbation travels across the bins.
 *              Shape:     Sine / Triangle / Pulse
 *              Direction: Forward / Backward / Alternating / Random
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
    float preCurve;          // 0..1 Hermite blend on base waveform before modulation
    float postSmooth;        // 0..1 global Hermite blend on final output
    bool  stepped;

    float phase;
    float phaseInc;
    float lastCvVal;     // cached CV value — only recompute phaseInc when CV changes

    PeriodicHermite16 herm;      // post-smooth hermite (on S[])
    PeriodicHermite16 preHerm;   // pre-smooth hermite (on baseRaw[])
    float baseRaw[kBins];        // raw base bins (morph only, no mod) for pre-smooth eval
    bool  preHermDirty;          // true when baseRaw has changed and slopes need recomputing
    bool  postHermDirty;         // true when S[] has changed and slopes need recomputing
    bool  modXDirty;             // true when kink/scramble has advanced and modX needs recomputing
    bool  sDirty;                // true when modX or bins have changed and S[] needs recomputing

    // Per-mode post-smooth
    bool  kinkStepped;      // stepped readout for kink shape only

    // Kink glide: slews kink position between steps
    float kinkGlide;             // 0..1 slew coefficient (0=instant, 1=never arrives)
    float kinkGlideSlewCoeff;    // cached per-sample slew
    float kinkPosSmoothed;       // current smoothed position

    // --- Modulation mode ---
    // 0=Off  1=Ripple  2=Wave  3=Scramble
    int   modMode;

    // --- Ripple mod ---
    float rippleAmt;
    float rippleLfoPhase;
    float rippleLfoPhaseInc;

    // --- Wave mod ---
    float kinkAmt;       // 0..1 depth
    int   kinkShape;     // 0=Triangle 1=Pulse
    int   kinkWidth;     // 1..16 bins, width of perturbation
    int   kinkRange;     // 1..16 bins, travel distance
    int   kinkDirection; // 0=Fwd 1=Back 2=Alt 3=Random
    int   kinkSampleCount;
    int   kinkSamplePeriod;

    // Wave position state — Kink 1
    float kinkPos;       // current centre position (float, 0..kBins-1)
    int   kinkAltDir;    // +1 or -1 for Alternating direction
    int   kinkRangeStart;// leftmost bin of travel range (set on first advance)
    bool  kinkRangeInit; // false until first advance initialises range

    // Kink 2 — shares Shape, Width, Range, Direction, Glide with K1
    float kink2Amt;
    int   kink2SampleCount;
    int   kink2SamplePeriod;
    float kink2Pos;
    float kink2PosSmoothed;
    float kink2GlideSlewCoeff;
    int   kink2AltDir;

    float modX[kBins];

    // --- Scramble mod ---
    float scrambleAmt;
    float scrambleGlide;         // 0..1 slew on modX between scramble steps
    float scrambleGlideSlewCoeff;// cached per-sample slew
    float scrambleModXTarget[kBins]; // target modX values (updated each step)
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

    // Recompute glide slew coefficient — call when Rate or Glide changes
    void updateKinkGlideSlewCoeff() {
        if (kinkGlide < 0.0001f) {
            kinkGlideSlewCoeff  = 1.0f;
            kink2GlideSlewCoeff = 1.0f;
        } else {
            float glideTime = kinkGlide * (float)kinkSamplePeriod;
            kinkGlideSlewCoeff = 1.0f - powf(0.001f, 1.0f / (glideTime + 1.0f));
            // K2 uses its own period for glide scaling
            float glideTime2 = kinkGlide * (float)kink2SamplePeriod;
            kink2GlideSlewCoeff = 1.0f - powf(0.001f, 1.0f / (glideTime2 + 1.0f));
        }
    }

    void updateScrambleGlideSlewCoeff() {
        if (scrambleGlide < 0.0001f) {
            scrambleGlideSlewCoeff = 1.0f;
        } else {
            float glideTime = scrambleGlide * (float)scrambleSamplePeriod;
            scrambleGlideSlewCoeff = 1.0f - powf(0.001f, 1.0f / (glideTime + 1.0f));
        }
    }
    // 1 → ~2s, 32 → ~25ms (exponential)
    static int rateToPeriod(int param, float sr) {
        float pct  = (float)param / 32.0f;
        float secs = 2.0f * powf(0.0125f / 2.0f, pct);
        int   p    = (int)(sr * secs);
        return (p < 1) ? 1 : p;
    }

    _CessnaAwg_DTC() {
        sampleRate = kDefaultSampleRate;
        depthD    = 0.5f;
        preCurve   = 0.0f;
        postSmooth = 0.0f;
        stepped    = false;

        kinkStepped       = false;
        kinkGlide         = 0.0f;
        kinkGlideSlewCoeff = 1.0f; // instant by default
        kinkPosSmoothed   = 0.0f;
        phase    = 0.0f;
        phaseInc = kC4Hz / kDefaultSampleRate;
        lastCvVal = -999.0f; // force recompute on first sample

        modMode = 0;

        // Ripple defaults
        rippleAmt          = 0.5f;
        rippleLfoPhase     = 0.0f;
        rippleLfoPhaseInc  = 0.5f / kDefaultSampleRate; // 0.5Hz default

        // Wave defaults
        kinkAmt         = 0.5f;
        kinkShape       = 0;    // Triangle
        kinkWidth       = 2;    // 2 bins wide
        kinkRange       = kBins; // full range
        kinkDirection   = 0;    // Forward
        kinkSampleCount = 0;
        kinkSamplePeriod = rateToPeriod(16, kDefaultSampleRate);
        kinkPos         = 0.0f;
        kinkAltDir      = 1;
        kinkRangeStart  = 0;
        kinkRangeInit   = false;

        // Kink 2 defaults — starts at opposite end of range
        kink2Amt              = 0.0f;  // off by default
        kink2SampleCount      = 0;
        kink2SamplePeriod     = rateToPeriod(12, kDefaultSampleRate); // slightly different default rate
        kink2Pos              = (float)(kBins - 1);
        kink2PosSmoothed      = (float)(kBins - 1);
        kink2GlideSlewCoeff   = 1.0f;
        kink2AltDir           = -1;  // starts going backward

        // Scramble defaults
        scrambleAmt          = 0.5f;
        scrambleGlide        = 0.0f;
        scrambleGlideSlewCoeff = 1.0f;
        scrambleSampleCount  = 0;
        scrambleSamplePeriod = rateToPeriod(16, kDefaultSampleRate);
        scrambleRng          = 0xABCD1234u;
        for (int i=0;i<kBins;i++) scrambleModXTarget[i] = 0.0f;

        for (int i=0;i<kBins;i++) {
            B[i] = S[i] = modX[i] = 0.0f;
            slotA[i] = slotB[i] = 0.0f;
            baseRaw[i] = 0.0f;
        }
        preHermDirty  = true;
        postHermDirty = true;
        modXDirty     = true;
        sDirty        = true;

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
    kParamPreCurve,     // pre-modulation interpolation blend
    kParamPostSmooth,   // post-modulation interpolation blend (global)
    kParamStepped,

    // Waveform page (bins 1-16)
    kParamB1,
    kParamB16 = kParamB1 + 15,

    // Mod selector page
    kParamModMode,

    // Ripple page
    kParamRippleAmt,
    kParamRippleRate,

    // Kink page
    kParamKinkAmt,       // K1 Amt
    kParamKinkRate,      // K1 Rate
    kParamKinkShape,
    kParamKinkWidth,
    kParamKinkRange,
    kParamKinkDirection,
    kParamKinkGlide,
    kParamKinkStepped,
    kParamKink2Amt,      // K2 Amt
    kParamKink2Rate,     // K2 Rate

    // Scramble page
    kParamScrambleAmt,
    kParamScrambleRate,
    kParamScrambleGlide,

    // Morph page
    kParamSaveA,
    kParamSaveB,
    kParamMorph,
};
static constexpr int kNumParams = kParamMorph + 1;

// ----------------------------
// String tables
// ----------------------------
static const char* modModeStrings[]  = {"Off","Ripple","Kink","Scramble",nullptr};
static const char* kinkShapeStrings[]= {"Triangle","Pulse",nullptr};
static const char* kinkDirStrings[]  = {"Forward","Backward","Alternating","Random",nullptr};
static const char* onOffStrings[]    = {"Off","On",nullptr};

// ----------------------------
// Parameter descriptors
// (must be in enum index order)
// ----------------------------
static const _NT_parameter g_parameters[kNumParams] = {
    NT_PARAMETER_AUDIO_INPUT("CV in", 0, 1)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("out", 0, 13)

    // Main
    { .name="Level",      .min=0,  .max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Pre-smooth", .min=0,  .max=100,.def=0,  .unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Post-smooth",.min=0,  .max=100,.def=100,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Stepped",    .min=0,  .max=1,  .def=0,  .unit=kNT_unitEnum,  .scaling=kNT_scalingNone,.enumStrings=onOffStrings },

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

    // Kink
    { .name="K1 Amt",    .min=0,.max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="K1 Rate",   .min=1,.max=32, .def=16,.unit=kNT_unitNone,  .scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Shape",     .min=0,.max=1,  .def=0, .unit=kNT_unitEnum,  .scaling=kNT_scalingNone,.enumStrings=kinkShapeStrings },
    { .name="Width",     .min=1,.max=16, .def=2, .unit=kNT_unitNone,  .scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Range",     .min=1,.max=16, .def=16,.unit=kNT_unitNone,  .scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Direction", .min=0,.max=3,  .def=0, .unit=kNT_unitEnum,  .scaling=kNT_scalingNone,.enumStrings=kinkDirStrings },
    { .name="Slew",      .min=0,.max=100,.def=0, .unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Angular",   .min=0,.max=1,  .def=0, .unit=kNT_unitEnum,  .scaling=kNT_scalingNone,.enumStrings=onOffStrings },
    { .name="K2 Amt",    .min=0,.max=100,.def=0, .unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="K2 Rate",   .min=1,.max=32, .def=12,.unit=kNT_unitNone,  .scaling=kNT_scalingNone,.enumStrings=nullptr },

    // Scramble
    { .name="Amt",        .min=0,.max=100,.def=50,.unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Rate",       .min=1,.max=32, .def=16,.unit=kNT_unitNone,  .scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Slew",       .min=0,.max=100,.def=0, .unit=kNT_unitPercent,.scaling=kNT_scalingNone,.enumStrings=nullptr },

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
    (uint8_t)kParamDepth,(uint8_t)kParamPreCurve,(uint8_t)kParamPostSmooth,(uint8_t)kParamStepped,
};
static const uint8_t g_pageBinsIdx[] = {
    (uint8_t)(kParamB1+0),(uint8_t)(kParamB1+1),(uint8_t)(kParamB1+2),(uint8_t)(kParamB1+3),
    (uint8_t)(kParamB1+4),(uint8_t)(kParamB1+5),(uint8_t)(kParamB1+6),(uint8_t)(kParamB1+7),
    (uint8_t)(kParamB1+8),(uint8_t)(kParamB1+9),(uint8_t)(kParamB1+10),(uint8_t)(kParamB1+11),
    (uint8_t)(kParamB1+12),(uint8_t)(kParamB1+13),(uint8_t)(kParamB1+14),(uint8_t)(kParamB1+15),
};
static const uint8_t g_pageModIdx[]      = { (uint8_t)kParamModMode };
static const uint8_t g_pageRippleIdx[]   = { (uint8_t)kParamRippleAmt, (uint8_t)kParamRippleRate };
static const uint8_t g_pageKinkIdx[]     = {
    (uint8_t)kParamKinkAmt,(uint8_t)kParamKinkRate,(uint8_t)kParamKinkShape,
    (uint8_t)kParamKinkWidth,(uint8_t)kParamKinkRange,(uint8_t)kParamKinkDirection,
    (uint8_t)kParamKinkGlide,(uint8_t)kParamKinkStepped,
    (uint8_t)kParamKink2Amt,(uint8_t)kParamKink2Rate,
};
static const uint8_t g_pageScrambleIdx[] = {
    (uint8_t)kParamScrambleAmt,(uint8_t)kParamScrambleRate,(uint8_t)kParamScrambleGlide,
};
static const uint8_t g_pageMorphIdx[]    = { (uint8_t)kParamSaveA,(uint8_t)kParamSaveB,(uint8_t)kParamMorph };

static const _NT_parameterPage g_pages_arr[] = {
    { .name="Routing", .numParams=(uint8_t)ARRAY_SIZE(g_pageRoutingIdx),.group=0,.unused={0,0},.params=g_pageRoutingIdx },
    { .name="Main",    .numParams=(uint8_t)ARRAY_SIZE(g_pageMainIdx),   .group=0,.unused={0,0},.params=g_pageMainIdx },
    { .name="Waveform",.numParams=(uint8_t)ARRAY_SIZE(g_pageBinsIdx),   .group=1,.unused={0,0},.params=g_pageBinsIdx },
    { .name="Mod",     .numParams=(uint8_t)ARRAY_SIZE(g_pageModIdx),    .group=0,.unused={0,0},.params=g_pageModIdx },
    { .name="Ripple",  .numParams=(uint8_t)ARRAY_SIZE(g_pageRippleIdx), .group=0,.unused={0,0},.params=g_pageRippleIdx },
    { .name="Kink",    .numParams=(uint8_t)ARRAY_SIZE(g_pageKinkIdx),   .group=0,.unused={0,0},.params=g_pageKinkIdx },
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
// combineS: slots + morph + pre-smooth + modX → S[]
// Pre-smooth: Hermite interpolation of base bins evaluated at each bin's
// centre phase position, blended with linear. Modulation then operates on
// this smoothed base rather than the raw stepped bin values.
// ----------------------------
// combineS: slots + morph + modX → S[]
// Pre-smooth is applied at audio rate in the step loop, not here.
// combineS just combines base + modX into S[] for the kink/mod to read.
// ----------------------------
static void combineS(struct _CessnaAwg_DTC *d) {
    if (!d->saveA) for (int i=0;i<kBins;i++) d->slotA[i] = d->B[i];
    if (!d->saveB) for (int i=0;i<kBins;i++) d->slotB[i] = d->B[i];

    float D = d->depthD;
    for (int i=0;i<kBins;i++) {
        float base = lerpf(d->slotA[i], d->slotB[i], d->morph);
        float v = base + d->modX[i] * 5.0f;
        float s = (v > 5.0f || v < -5.0f) ? tanhf(v / 5.0f) * 5.0f : v;
        d->S[i] = D * s;
    }

    // Store raw base for pre-smooth evaluation at audio rate
    for (int i=0;i<kBins;i++)
        d->baseRaw[i] = lerpf(d->slotA[i], d->slotB[i], d->morph);
    d->preHermDirty  = true;
    d->postHermDirty = true;
}

// ----------------------------
// Wave modulation
// ----------------------------

// Build the perturbation shape envelope for a given position
// Returns values in 0..1 (unsigned) — polarity applied by caller
static float kinkEnvelope(int shape, int binIdx, int centrePos, int halfWidth) {
    int dist = binIdx - centrePos;
    // Wrap distance to nearest (handles wraparound)
    if (dist >  kBins/2) dist -= kBins;
    if (dist < -kBins/2) dist += kBins;
    int absDist = dist < 0 ? -dist : dist;

    if (absDist > halfWidth) return 0.0f;

    // Guard against division by zero when halfWidth=0 (Width=1)
    // — both shapes collapse to a single bin spike
    if (halfWidth == 0) return (absDist == 0) ? 1.0f : 0.0f;

    float t = (float)absDist / (float)halfWidth;

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
static void advanceKinkPosition(struct _CessnaAwg_DTC *d) {
    // Initialise range anchor on first advance
    if (!d->kinkRangeInit) {
        d->kinkRangeStart = 0;
        d->kinkPos        = 0.0f;
        d->kinkRangeInit  = true;
    }

    float range = (float)(d->kinkRange - 1);
    if (range < 0.0f) range = 0.0f;

    switch (d->kinkDirection) {
        case 0: // Forward — step +1, wrap within range
            d->kinkPos = d->kinkPos + 1.0f;
            if (d->kinkPos > range) d->kinkPos = 0.0f;
            break;
        case 1: // Backward — step -1, wrap within range
            d->kinkPos = d->kinkPos - 1.0f;
            if (d->kinkPos < 0.0f) d->kinkPos = range;
            break;
        case 2: // Alternating — pingpong within range
            d->kinkPos += (float)d->kinkAltDir;
            if (d->kinkPos >= range) { d->kinkPos = range; d->kinkAltDir = -1; }
            if (d->kinkPos <= 0.0f) { d->kinkPos = 0.0f;  d->kinkAltDir =  1; }
            break;
        case 3: { // Random — jump to random position within range
            float r = lcgRand01(d->scrambleRng);
            d->kinkPos = r * range;
            break;
        }
    }
}

// Advance Kink 2 position — shares direction logic with K1 but independent state
static void advanceKink2Position(struct _CessnaAwg_DTC *d) {
    float range = (float)(d->kinkRange - 1);
    if (range < 0.0f) range = 0.0f;

    switch (d->kinkDirection) {
        case 0: // Forward
            d->kink2Pos = d->kink2Pos + 1.0f;
            if (d->kink2Pos > range) d->kink2Pos = 0.0f;
            break;
        case 1: // Backward
            d->kink2Pos = d->kink2Pos - 1.0f;
            if (d->kink2Pos < 0.0f) d->kink2Pos = range;
            break;
        case 2: // Alternating — starts in opposite direction to K1
            d->kink2Pos += (float)d->kink2AltDir;
            if (d->kink2Pos >= range) { d->kink2Pos = range; d->kink2AltDir = -1; }
            if (d->kink2Pos <= 0.0f) { d->kink2Pos = 0.0f;  d->kink2AltDir =  1; }
            break;
        case 3: { // Random
            float r = lcgRand01(d->scrambleRng);
            d->kink2Pos = r * range;
            break;
        }
    }
}

// Update kink position slew every sample — returns true if centrePos changed
static bool updateKinkSlew(struct _CessnaAwg_DTC *d) {
    float target1 = (float)(((int)roundf(d->kinkPos) + d->kinkRangeStart + kBins * 4) % kBins);
    float target2 = (float)(((int)roundf(d->kink2Pos) + d->kinkRangeStart + kBins * 4) % kBins);

    int oldCentre1 = ((int)roundf(d->kinkPosSmoothed) % kBins + kBins) % kBins;
    int oldCentre2 = ((int)roundf(d->kink2PosSmoothed) % kBins + kBins) % kBins;

    if (d->kinkGlide > 0.0001f) {
        float diff1 = target1 - d->kinkPosSmoothed;
        if (diff1 >  (float)(kBins/2)) diff1 -= (float)kBins;
        if (diff1 < -(float)(kBins/2)) diff1 += (float)kBins;
        d->kinkPosSmoothed += d->kinkGlideSlewCoeff * diff1;
        if (d->kinkPosSmoothed < 0.0f)         d->kinkPosSmoothed += (float)kBins;
        if (d->kinkPosSmoothed >= (float)kBins) d->kinkPosSmoothed -= (float)kBins;

        if (d->kink2Amt > 0.0001f) {
            float diff2 = target2 - d->kink2PosSmoothed;
            if (diff2 >  (float)(kBins/2)) diff2 -= (float)kBins;
            if (diff2 < -(float)(kBins/2)) diff2 += (float)kBins;
            d->kink2PosSmoothed += d->kink2GlideSlewCoeff * diff2;
            if (d->kink2PosSmoothed < 0.0f)         d->kink2PosSmoothed += (float)kBins;
            if (d->kink2PosSmoothed >= (float)kBins) d->kink2PosSmoothed -= (float)kBins;
        }
    } else {
        d->kinkPosSmoothed  = target1;
        d->kink2PosSmoothed = target2;
    }

    int newCentre1 = ((int)roundf(d->kinkPosSmoothed) % kBins + kBins) % kBins;
    int newCentre2 = ((int)roundf(d->kink2PosSmoothed) % kBins + kBins) % kBins;

    return (newCentre1 != oldCentre1) || (newCentre2 != oldCentre2);
}

// Build modX from current smoothed kink positions — sums K1 and K2 contributions
static void buildKinkModX(struct _CessnaAwg_DTC *d) {
    int centrePos1 = ((int)roundf(d->kinkPosSmoothed) % kBins + kBins) % kBins;
    int centrePos2 = ((int)roundf(d->kink2PosSmoothed) % kBins + kBins) % kBins;

    int halfWidth = d->kinkWidth / 2;
    if (halfWidth < 1) halfWidth = 1;
    int halfWidth1 = (d->kinkShape == 1) ? (halfWidth / 2) : halfWidth;
    if (halfWidth1 < 1) halfWidth1 = 1;

    for (int i=0;i<kBins;i++) {
        float env1 = kinkEnvelope(d->kinkShape, i, centrePos1, halfWidth1);
        float env2 = (d->kink2Amt > 0.0001f)
                   ? kinkEnvelope(d->kinkShape, i, centrePos2, halfWidth1)
                   : 0.0f;
        d->modX[i] = d->kinkAmt * env1 + d->kink2Amt * env2;
    }
}

// ----------------------------
// Scramble modulation
// ----------------------------
static void advanceScramble(struct _CessnaAwg_DTC *d) {
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
    // Store as target — modX will slew toward this
    for (int i=0;i<kBins;i++)
        d->scrambleModXTarget[i] = (tmp[i] - d->B[i]) / 5.0f;
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
        case 2: // Wave — modX already built by advanceKinkPosition/buildKinkModX
            buildKinkModX(d);
            break;
        case 3: // Scramble — slew modX toward scrambleModXTarget
            for (int i=0;i<kBins;i++)
                d->modX[i] += d->scrambleGlideSlewCoeff * (d->scrambleModXTarget[i] - d->modX[i]);
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
        case kParamDepth:
            d->depthD = getPct01(kParamDepth);
            d->sDirty = true;
            break;
        case kParamPreCurve:   d->preCurve   = getPct01(kParamPreCurve); break;
        case kParamPostSmooth: d->postSmooth = getPct01(kParamPostSmooth); break;
        case kParamStepped:    d->stepped    = (self->v[kParamStepped] != 0); break;

        // Mod selector
        case kParamModMode:
            d->modMode = self->v[kParamModMode];
            if (d->modMode == 2) { // Kink
                d->kinkPos          = 0.0f;
                d->kinkPosSmoothed  = 0.0f;
                d->kinkRangeInit    = false;
                d->kinkAltDir       = 1;
                float range2 = (float)(d->kinkRange - 1);
                d->kink2Pos         = range2;
                d->kink2PosSmoothed = range2;
                d->kink2AltDir      = -1;
                d->kink2SampleCount = 0;
            }
            d->modXDirty = true;
            d->sDirty    = true;
            break;

        // Ripple
        case kParamRippleAmt:
            d->rippleAmt = getPct01(kParamRippleAmt);
            d->modXDirty = true;
            break;
        case kParamRippleRate: {
            float hz = (float)self->v[kParamRippleRate] * 0.01f;
            d->rippleLfoPhaseInc = hz / d->sampleRate;
            break;
        }

        // Wave
        case kParamKinkAmt:
            d->kinkAmt = getPct01(kParamKinkAmt);
            d->modXDirty = true;
            break;
        case kParamKinkRate:
            d->kinkSamplePeriod = _CessnaAwg_DTC::rateToPeriod(self->v[kParamKinkRate], d->sampleRate);
            d->updateKinkGlideSlewCoeff();
            break;
        case kParamKinkShape:     d->kinkShape     = self->v[kParamKinkShape]; break;
        case kParamKinkWidth:     d->kinkWidth     = self->v[kParamKinkWidth]; break;
        case kParamKinkRange:
            d->kinkRange     = self->v[kParamKinkRange];
            if (d->kinkPos >= (float)d->kinkRange) d->kinkPos = 0.0f;
            break;
        case kParamKinkDirection:
            d->kinkDirection = self->v[kParamKinkDirection];
            break;
        case kParamKinkGlide:
            d->kinkGlide = getPct01(kParamKinkGlide);
            d->updateKinkGlideSlewCoeff();
            break;
            break;
        case kParamKinkStepped:
            d->kinkStepped = (self->v[kParamKinkStepped] != 0);
            break;
        case kParamKink2Amt:
            d->kink2Amt = getPct01(kParamKink2Amt);
            d->modXDirty = true;
            break;
        case kParamKink2Rate:
            d->kink2SamplePeriod = _CessnaAwg_DTC::rateToPeriod(self->v[kParamKink2Rate], d->sampleRate);
            d->updateKinkGlideSlewCoeff();
            break;

        // Scramble
        case kParamScrambleAmt:
            d->scrambleAmt = getPct01(kParamScrambleAmt);
            d->modXDirty = true;
            break;
        case kParamScrambleRate:
            d->scrambleSamplePeriod = _CessnaAwg_DTC::rateToPeriod(self->v[kParamScrambleRate], d->sampleRate);
            d->updateScrambleGlideSlewCoeff();
            break;
        case kParamScrambleGlide:
            d->scrambleGlide = getPct01(kParamScrambleGlide);
            d->updateScrambleGlideSlewCoeff();
            break;
            break;

        // Morph
        case kParamSaveA: {
            d->saveA = (self->v[kParamSaveA] != 0);
            if (d->saveA) for (int i=0;i<kBins;i++) d->slotA[i] = d->B[i];
            d->sDirty = true;
            d->preHermDirty = true;
            break;
        }
        case kParamSaveB: {
            d->saveB = (self->v[kParamSaveB] != 0);
            if (d->saveB) for (int i=0;i<kBins;i++) d->slotB[i] = d->B[i];
            d->sDirty = true;
            d->preHermDirty = true;
            break;
        }
        case kParamMorph:
            d->morph = getPct01(kParamMorph);
            d->sDirty = true;
            d->preHermDirty = true;
            break;

        default:
            if (p >= kParamB1 && p <= kParamB16) {
                d->B[p - kParamB1] = getPctSigned(p);
                d->sDirty      = true;
                d->preHermDirty = true;
            }
            break;
    }
}

// ----------------------------
// updateDisplay
// ----------------------------
static void updateDisplay(struct _CessnaAwg_DTC *d) {
    float postCurve = d->postSmooth;

    bool useGlobalStep = d->stepped;
    bool useKinkStep   = (d->modMode == 2 && d->kinkStepped);

    if (!useGlobalStep && d->preCurve > 0.0001f && d->preHermDirty) {
        d->preHerm.computeSlopes(d->baseRaw);
        d->preHermDirty = false;
    }
    if (!useGlobalStep && !useKinkStep && postCurve > 0.0001f && d->postHermDirty) {
        d->herm.computeSlopes(d->S);
        d->postHermDirty = false;
    }

    for (int k=0;k<kDisplayPoints;k++) {
        float ph = (float)k / (float)(kDisplayPoints-1);
        float y;
        if (useGlobalStep) {
            y = evalStepped(d->S, ph);
        } else if (useKinkStep) {
            float baseLin  = evalLinear(d->baseRaw, ph);
            float baseHerm = (d->preCurve > 0.0001f) ? d->preHerm.eval(d->baseRaw, ph) : baseLin;
            float baseVal  = lerpf(baseLin, baseHerm, d->preCurve) * d->depthD;
            float kinkVal  = evalStepped(d->modX, ph) * 5.0f * d->depthD;
            float combined = baseVal + kinkVal;
            y = (combined > 5.0f || combined < -5.0f) ? tanhf(combined / 5.0f) * 5.0f : combined;
        } else if (postCurve > 0.0001f) {
            float sLin  = evalLinear(d->S, ph);
            float sHerm = d->herm.eval(d->S, ph);
            y = lerpf(sLin, sHerm, postCurve);
        } else {
            float baseLin  = evalLinear(d->baseRaw, ph);
            float baseHerm = (d->preCurve > 0.0001f) ? d->preHerm.eval(d->baseRaw, ph) : baseLin;
            float baseVal  = lerpf(baseLin, baseHerm, d->preCurve) * d->depthD;
            float modVal   = evalLinear(d->modX, ph) * 5.0f * d->depthD;
            float combined = baseVal + modVal;
            y = (combined > 5.0f || combined < -5.0f) ? tanhf(combined / 5.0f) * 5.0f : combined;
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
        // Pitch — only recompute powf when CV changes meaningfully
        float cvVal = cv ? cv[n] : 0.0f;
        if (fabsf(cvVal - d->lastCvVal) > 0.0001f) {
            d->lastCvVal = cvVal;
            d->phaseInc  = kC4Hz * powf(2.0f, cvVal) / d->sampleRate;
        }

        // Timed advances for Kink
        if (d->modMode == 2) {
            // Kink 1 advance
            d->kinkSampleCount++;
            if (d->kinkSampleCount >= d->kinkSamplePeriod) {
                d->kinkSampleCount = 0;
                advanceKinkPosition(d);
            }
            // Kink 2 advance (only if active)
            if (d->kink2Amt > 0.0001f) {
                d->kink2SampleCount++;
                if (d->kink2SampleCount >= d->kink2SamplePeriod) {
                    d->kink2SampleCount = 0;
                    advanceKink2Position(d);
                }
            }
            // Slew runs every sample — only dirty if centrePos changed
            if (updateKinkSlew(d))
                d->modXDirty = true;
        } else if (d->modMode == 3) {
            // Scramble: timed advance
            d->scrambleSampleCount++;
            if (d->scrambleSampleCount >= d->scrambleSamplePeriod) {
                d->scrambleSampleCount = 0;
                advanceScramble(d);
            }
            // Scramble slew is continuous — always dirty
            d->modXDirty = true;
        } else if (d->modMode == 1) {
            // Ripple is continuous — always dirty
            d->modXDirty = true;
        }

        // Compute modX only when needed
        if (d->modXDirty) {
            computeModX(d);
            d->modXDirty = false;
            d->sDirty    = true;
        }

        // Combine → S[] only when modX or bins have changed
        if (d->sDirty) {
            combineS(d);
            d->sDirty        = false;
            d->postHermDirty = true;
        }

        // Post-smooth: global + per-mode (take max)
        float postCurve = d->postSmooth;

        // Determine if fully stepped (global) or kink-only stepped
        bool useGlobalStep  = d->stepped;
        bool useKinkStep    = (d->modMode == 2 && d->kinkStepped);

        // Pre-smooth slopes — only recompute when baseRaw has changed
        if (!useGlobalStep && d->preCurve > 0.0001f && d->preHermDirty) {
            d->preHerm.computeSlopes(d->baseRaw);
            d->preHermDirty = false;
        }
        // Post-smooth slopes — only recompute when S[] has changed
        if (!useGlobalStep && !useKinkStep && postCurve > 0.0001f && d->postHermDirty) {
            d->herm.computeSlopes(d->S);
            d->postHermDirty = false;
        }

        // Evaluate wavetable
        float y;
        if (useGlobalStep) {
            // Whole waveform stepped
            y = evalStepped(d->S, d->phase);
        } else if (useKinkStep) {
            // Base smooth, kink stepped
            float baseLin  = evalLinear(d->baseRaw, d->phase);
            float baseHerm = (d->preCurve > 0.0001f) ? d->preHerm.eval(d->baseRaw, d->phase) : baseLin;
            float baseVal  = lerpf(baseLin, baseHerm, d->preCurve) * d->depthD;
            float kinkVal  = evalStepped(d->modX, d->phase) * 5.0f * d->depthD;
            float combined = baseVal + kinkVal;
            y = (combined > 5.0f || combined < -5.0f) ? tanhf(combined / 5.0f) * 5.0f : combined;
        } else {
            // Base: blend linear and Hermite using preCurve
            float baseLin  = evalLinear(d->baseRaw, d->phase);
            float baseHerm = (d->preCurve > 0.0001f) ? d->preHerm.eval(d->baseRaw, d->phase) : baseLin;
            float baseVal  = lerpf(baseLin, baseHerm, d->preCurve) * d->depthD;
            float modVal   = evalLinear(d->modX, d->phase) * 5.0f * d->depthD;
            float combined = baseVal + modVal;
            if (combined > 5.0f || combined < -5.0f)
                combined = tanhf(combined / 5.0f) * 5.0f;
            // Post-smooth
            if (postCurve > 0.0001f) {
                float sLin  = evalLinear(d->S, d->phase);
                float sHerm = d->herm.eval(d->S, d->phase);
                y = lerpf(sLin, sHerm, postCurve);
            } else {
                y = combined;
            }
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
    d->kinkSamplePeriod      = _CessnaAwg_DTC::rateToPeriod(16, sr);
    d->kink2SamplePeriod     = _CessnaAwg_DTC::rateToPeriod(12, sr);
    d->updateKinkGlideSlewCoeff();
    d->scrambleSamplePeriod  = _CessnaAwg_DTC::rateToPeriod(16, sr);
    d->updateScrambleGlideSlewCoeff();

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
