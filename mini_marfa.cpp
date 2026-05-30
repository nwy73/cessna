/*
 * MiniMARFA v1 — 8-stage MARF-style arbitrary function programmer
 * -----------------------------------------------------------------
 * A single 8-stage AFG, reduced from the Buchla 248/MARF concept for
 * disting NT.  The intention is not to modernise the MARF into a normal
 * clocked sequencer, but to preserve the MARF's basic semantics within
 * NT constraints.
 *
 * Core semantics from the MARF manual:
 *   - Pulse outputs fire when a stage is addressed / at stage entry.
 *   - All Pulses fires at the start of every stage.
 *   - STOP stops the internal clock once it reaches that stage.
 *   - SUST holds when the stage is reached as long as START is high.
 *   - ENABLE holds when the stage is reached unless START is high.
 *   - FIRST/LAST define cycles.
 *   - No FIRST/LAST markers: free-run through all 8 stages.
 *   - QUANT samples at stage entry and quantizes to 0.1V increments.
 *   - CONT follows live CV level changes.
 *   - SLOPED ramps to the stage value.
 *   - STEPPED jumps immediately to the stage value.
 *
 * v1 simplifications:
 *   - One AFG only.
 *   - 8 stages only.
 *   - Quant/Cont, Sloped/Stepped, Time Range, and Voltage Range are global.
 *   - Stage addressing is included, but simple: STROB + SEXT jumps to a stage;
 *     continuous external addressing immediately follows SEXT.
 */

#include <math.h>
#include <string.h>
#include <new>
#include <distingnt/api.h>

// Gallery/public UUID placeholder — replace before publishing if desired.
// Plugin GUID: 2D7699AB-3C80-4D26-9C65-24EED393C4A8

extern "C" int* __errno(void) {
    static int errno_val = 0;
    return &errno_val;
}

#ifndef ARRAY_SIZE
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))
#endif

// ----------------------------
// Constants
// ----------------------------
static constexpr int   kStages            = 8;
static constexpr int   kDisplayPoints     = 128;
static constexpr int   kMaxBus            = 28;
static constexpr float kDefaultSampleRate = 48000.0f;
static constexpr float kGateThreshold     = 1.0f;
static constexpr float kPulseVolts        = 10.0f;

// ----------------------------
// Helpers
// ----------------------------
static inline float clampf(float x, float lo, float hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}
static inline float lerpf(float a, float b, float t) {
    return a + (b - a) * t;
}
static inline bool risingEdge(float v, bool &oldHigh) {
    bool high = (v > kGateThreshold);
    bool edge = high && !oldHigh;
    oldHigh = high;
    return edge;
}

// ----------------------------
// Stage flags — order follows MARF panel/programming order
// ----------------------------
struct StageFlags {
    bool pulse1;
    bool pulse2;
    bool stop;
    bool sust;
    bool enable;
    bool first;
    bool last;

    StageFlags()
    : pulse1(false), pulse2(false), stop(false), sust(false),
      enable(false), first(false), last(false) {}
};

// ----------------------------
// DTC
// ----------------------------
struct _MiniMARFA_DTC {
    float sampleRate;

    StageFlags flags[kStages];

    int   currentStage;
    int   previousStage;
    int   cycleFirst;
    int   cycleLast;

    bool  running;
    bool  held;
    bool  startGate;
    bool  stopHigh;
    bool  resetHigh;
    bool  strobeHigh;

    float phase;                 // 0..1 through current stage
    float stageDurationSeconds;
    float output;
    float fromVoltage;
    float targetVoltage;
    float sampledVoltage;

    int   pulse1Remaining;
    int   pulse2Remaining;
    int   allPulseRemaining;

    float dispY[kDisplayPoints];
    int   dispStageX[kStages];

    _MiniMARFA_DTC() {
        sampleRate = kDefaultSampleRate;
        currentStage = 0;
        previousStage = 0;
        cycleFirst = 0;
        cycleLast = kStages - 1;

        running = false;
        held = true;
        startGate = stopHigh = resetHigh = strobeHigh = false;

        phase = 0.0f;
        stageDurationSeconds = 0.1f;
        output = 0.0f;
        fromVoltage = 0.0f;
        targetVoltage = 0.0f;
        sampledVoltage = 0.0f;

        pulse1Remaining = pulse2Remaining = allPulseRemaining = 0;

        for (int i=0;i<kDisplayPoints;i++) dispY[i] = 0.0f;
        for (int i=0;i<kStages;i++) dispStageX[i] = 0;
    }
};

// ----------------------------
// Plugin struct
// ----------------------------
struct _MiniMARFA : public _NT_algorithm {
    _MiniMARFA_DTC *dtc;
    _MiniMARFA(_MiniMARFA_DTC *d) : dtc(d) {}
};

// ----------------------------
// Parameter enum
// ----------------------------
enum {
    // Routing page
    kParamVOut = 0, kParamVOutMode,
    kParamTOut,     kParamTOutMode,
    kParamROut,     kParamROutMode,
    kParamP1Out,    kParamP1OutMode,
    kParamP2Out,    kParamP2OutMode,
    kParamAPOut,    kParamAPOutMode,
    kParamStartIn,
    kParamStopIn,
    kParamResetIn,
    kParamStrobeIn,
    kParamSExtIn,

    // CV Levels page
    kParamCV1,
    kParamCV8 = kParamCV1 + 7,

    // Time Levels page
    kParamTime1,
    kParamTime8 = kParamTime1 + 7,

    // Per-stage program pages (8 stages × 7 flags = 56 params)
    // All 56 params are contiguous: kParamS1Pulse1 + stage*kFlagsPerStage + kFlagXxx
    kParamS1Pulse1,
    kParamS8Last   = kParamS1Pulse1 + 8*7 - 1,

    // Global page — starts immediately after the 56 stage flag params
    kParamManualStart  = kParamS8Last + 1,
    kParamManualStop,
    kParamManualReset,
    kParamManualStrobe,
    kParamQuantCont,
    kParamSlopedStepped,
    kParamTimeRange,
    kParamVoltageRange,
    kParamStageAddress,
    kParamPulseLength,
};
static constexpr int kNumParams = kParamPulseLength + 1;
static constexpr int kFlagsPerStage = 7;  // Pulse1, Pulse2, Stop, Sust, Enable, First, Last
static constexpr int kFlagPulse1 = 0, kFlagPulse2 = 1, kFlagStop = 2,
                     kFlagSust  = 3, kFlagEnable = 4, kFlagFirst = 5, kFlagLast = 6;

// ----------------------------
// String tables
// ----------------------------
static const char* onOffStrings[]          = {"Off", "On", nullptr};
static const char* quantContStrings[]      = {"Quant", "Cont", nullptr};
static const char* slopedSteppedStrings[]  = {"Sloped", "Stepped", nullptr};
static const char* timeRangeStrings[]      = {".002-.03s", ".02-.3s", ".2-3s", "2-30s", nullptr};
static const char* voltageRangeStrings[]   = {"0-10V", "0-5V", "+/-5V", nullptr};
static const char* stageAddressStrings[]   = {"Internal", "Strobe Ext", "Cont Ext", nullptr};

// ----------------------------
// Parameter descriptors
// ----------------------------
static const _NT_parameter g_parameters[kNumParams] = {
    // Routing. Outputs use WITH_MODE to match the working Cessna pattern.
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Function CV", 1, 13)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Time",        0, 0)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Ref",         0, 0)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Pulse 1",     0, 0)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Pulse 2",     0, 0)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("All Pulses",  0, 0)
    NT_PARAMETER_AUDIO_INPUT("Start",  0, 1)
    NT_PARAMETER_AUDIO_INPUT("Stop",   0, 2)
    NT_PARAMETER_AUDIO_INPUT("Reset",  0, 0)
    NT_PARAMETER_AUDIO_INPUT("Strobe", 0, 0)
    NT_PARAMETER_AUDIO_INPUT("Stage Ext", 0, 0)

    // CV levels — 0..1000 maps to selected voltage range.
    { .name="CV 1", .min=0,.max=1000,.def=0,   .unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="CV 2", .min=0,.max=1000,.def=1000,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="CV 3", .min=0,.max=1000,.def=500, .unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="CV 4", .min=0,.max=1000,.def=750, .unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="CV 5", .min=0,.max=1000,.def=250, .unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="CV 6", .min=0,.max=1000,.def=1000,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="CV 7", .min=0,.max=1000,.def=500, .unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="CV 8", .min=0,.max=1000,.def=0,   .unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },

    // Time levels — 0..1000 within selected global range.
    { .name="Time 1", .min=0,.max=1000,.def=500,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Time 2", .min=0,.max=1000,.def=500,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Time 3", .min=0,.max=1000,.def=500,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Time 4", .min=0,.max=1000,.def=500,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Time 5", .min=0,.max=1000,.def=500,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Time 6", .min=0,.max=1000,.def=500,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Time 7", .min=0,.max=1000,.def=500,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Time 8", .min=0,.max=1000,.def=500,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },

    // Per-stage program params (8 stages × 7 flags = 56 params).
    // Indexed as kParamS1Pulse1 + stage*kFlagsPerStage + kFlagXxx.
#define STAGE_FLAGS \
    { .name="Pulse 1",.min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings }, \
    { .name="Pulse 2",.min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings }, \
    { .name="Stop",   .min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings }, \
    { .name="Sust",   .min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings }, \
    { .name="Enable", .min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings }, \
    { .name="First",  .min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings }, \
    { .name="Last",   .min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings },
    STAGE_FLAGS STAGE_FLAGS STAGE_FLAGS STAGE_FLAGS
    STAGE_FLAGS STAGE_FLAGS STAGE_FLAGS STAGE_FLAGS
#undef STAGE_FLAGS


    // Global — manual start/stop buttons, then mode settings
    { .name="Start",         .min=0,.max=1,.def=0,.unit=kNT_unitConfirm,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Stop",          .min=0,.max=1,.def=0,.unit=kNT_unitConfirm,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Reset",         .min=0,.max=1,.def=0,.unit=kNT_unitConfirm,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Strobe",        .min=0,.max=1,.def=0,.unit=kNT_unitConfirm,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Quant/Cont",    .min=0,.max=1,.def=1,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=quantContStrings },
    { .name="Sloped/Stepped",.min=0,.max=1,.def=1,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=slopedSteppedStrings },
    { .name="Time Range",    .min=0,.max=3,.def=1,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=timeRangeStrings },
    { .name="Voltage Range", .min=0,.max=2,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=voltageRangeStrings },
    { .name="Stage Address", .min=0,.max=2,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=stageAddressStrings },
    { .name="Pulse ms",      .min=1,.max=100,.def=1,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
};

// ----------------------------
// Pages
// ----------------------------
static const uint8_t g_pageRoutingIdx[] = {
    (uint8_t)kParamVOut,(uint8_t)kParamVOutMode,
    (uint8_t)kParamTOut,(uint8_t)kParamTOutMode,
    (uint8_t)kParamROut,(uint8_t)kParamROutMode,
    (uint8_t)kParamP1Out,(uint8_t)kParamP1OutMode,
    (uint8_t)kParamP2Out,(uint8_t)kParamP2OutMode,
    (uint8_t)kParamAPOut,(uint8_t)kParamAPOutMode,
    (uint8_t)kParamStartIn,(uint8_t)kParamStopIn,
    (uint8_t)kParamResetIn,(uint8_t)kParamStrobeIn,(uint8_t)kParamSExtIn,
};
static const uint8_t g_pageCVIdx[] = {
    (uint8_t)(kParamCV1+0),(uint8_t)(kParamCV1+1),(uint8_t)(kParamCV1+2),(uint8_t)(kParamCV1+3),
    (uint8_t)(kParamCV1+4),(uint8_t)(kParamCV1+5),(uint8_t)(kParamCV1+6),(uint8_t)(kParamCV1+7),
};
static const uint8_t g_pageTimeIdx[] = {
    (uint8_t)(kParamTime1+0),(uint8_t)(kParamTime1+1),(uint8_t)(kParamTime1+2),(uint8_t)(kParamTime1+3),
    (uint8_t)(kParamTime1+4),(uint8_t)(kParamTime1+5),(uint8_t)(kParamTime1+6),(uint8_t)(kParamTime1+7),
};
static const uint8_t g_pageStage1Idx[] = {
    (uint8_t)(kParamS1Pulse1+0),(uint8_t)(kParamS1Pulse1+1),(uint8_t)(kParamS1Pulse1+2),
    (uint8_t)(kParamS1Pulse1+3),(uint8_t)(kParamS1Pulse1+4),(uint8_t)(kParamS1Pulse1+5),
    (uint8_t)(kParamS1Pulse1+6),
};
static const uint8_t g_pageStage2Idx[] = {
    (uint8_t)(kParamS1Pulse1+7),(uint8_t)(kParamS1Pulse1+8),(uint8_t)(kParamS1Pulse1+9),
    (uint8_t)(kParamS1Pulse1+10),(uint8_t)(kParamS1Pulse1+11),(uint8_t)(kParamS1Pulse1+12),
    (uint8_t)(kParamS1Pulse1+13),
};
static const uint8_t g_pageStage3Idx[] = {
    (uint8_t)(kParamS1Pulse1+14),(uint8_t)(kParamS1Pulse1+15),(uint8_t)(kParamS1Pulse1+16),
    (uint8_t)(kParamS1Pulse1+17),(uint8_t)(kParamS1Pulse1+18),(uint8_t)(kParamS1Pulse1+19),
    (uint8_t)(kParamS1Pulse1+20),
};
static const uint8_t g_pageStage4Idx[] = {
    (uint8_t)(kParamS1Pulse1+21),(uint8_t)(kParamS1Pulse1+22),(uint8_t)(kParamS1Pulse1+23),
    (uint8_t)(kParamS1Pulse1+24),(uint8_t)(kParamS1Pulse1+25),(uint8_t)(kParamS1Pulse1+26),
    (uint8_t)(kParamS1Pulse1+27),
};
static const uint8_t g_pageStage5Idx[] = {
    (uint8_t)(kParamS1Pulse1+28),(uint8_t)(kParamS1Pulse1+29),(uint8_t)(kParamS1Pulse1+30),
    (uint8_t)(kParamS1Pulse1+31),(uint8_t)(kParamS1Pulse1+32),(uint8_t)(kParamS1Pulse1+33),
    (uint8_t)(kParamS1Pulse1+34),
};
static const uint8_t g_pageStage6Idx[] = {
    (uint8_t)(kParamS1Pulse1+35),(uint8_t)(kParamS1Pulse1+36),(uint8_t)(kParamS1Pulse1+37),
    (uint8_t)(kParamS1Pulse1+38),(uint8_t)(kParamS1Pulse1+39),(uint8_t)(kParamS1Pulse1+40),
    (uint8_t)(kParamS1Pulse1+41),
};
static const uint8_t g_pageStage7Idx[] = {
    (uint8_t)(kParamS1Pulse1+42),(uint8_t)(kParamS1Pulse1+43),(uint8_t)(kParamS1Pulse1+44),
    (uint8_t)(kParamS1Pulse1+45),(uint8_t)(kParamS1Pulse1+46),(uint8_t)(kParamS1Pulse1+47),
    (uint8_t)(kParamS1Pulse1+48),
};
static const uint8_t g_pageStage8Idx[] = {
    (uint8_t)(kParamS1Pulse1+49),(uint8_t)(kParamS1Pulse1+50),(uint8_t)(kParamS1Pulse1+51),
    (uint8_t)(kParamS1Pulse1+52),(uint8_t)(kParamS1Pulse1+53),(uint8_t)(kParamS1Pulse1+54),
    (uint8_t)(kParamS1Pulse1+55),
};
static const uint8_t g_pageGlobalIdx[] = {
    (uint8_t)kParamManualStart,(uint8_t)kParamManualStop,
    (uint8_t)kParamManualReset,(uint8_t)kParamManualStrobe,
    (uint8_t)kParamQuantCont,(uint8_t)kParamSlopedStepped,(uint8_t)kParamTimeRange,
    (uint8_t)kParamVoltageRange,(uint8_t)kParamStageAddress,(uint8_t)kParamPulseLength,
};

static const _NT_parameterPage g_pages_arr[] = {
    { .name="Routing",     .numParams=(uint8_t)ARRAY_SIZE(g_pageRoutingIdx),.group=0,.unused={0,0},.params=g_pageRoutingIdx },
    { .name="CV Levels",   .numParams=(uint8_t)ARRAY_SIZE(g_pageCVIdx),     .group=1,.unused={0,0},.params=g_pageCVIdx },
    { .name="Time Levels", .numParams=(uint8_t)ARRAY_SIZE(g_pageTimeIdx),   .group=1,.unused={0,0},.params=g_pageTimeIdx },
    { .name="Stage 1",     .numParams=(uint8_t)ARRAY_SIZE(g_pageStage1Idx), .group=2,.unused={0,0},.params=g_pageStage1Idx },
    { .name="Stage 2",     .numParams=(uint8_t)ARRAY_SIZE(g_pageStage2Idx), .group=2,.unused={0,0},.params=g_pageStage2Idx },
    { .name="Stage 3",     .numParams=(uint8_t)ARRAY_SIZE(g_pageStage3Idx), .group=2,.unused={0,0},.params=g_pageStage3Idx },
    { .name="Stage 4",     .numParams=(uint8_t)ARRAY_SIZE(g_pageStage4Idx), .group=2,.unused={0,0},.params=g_pageStage4Idx },
    { .name="Stage 5",     .numParams=(uint8_t)ARRAY_SIZE(g_pageStage5Idx), .group=2,.unused={0,0},.params=g_pageStage5Idx },
    { .name="Stage 6",     .numParams=(uint8_t)ARRAY_SIZE(g_pageStage6Idx), .group=2,.unused={0,0},.params=g_pageStage6Idx },
    { .name="Stage 7",     .numParams=(uint8_t)ARRAY_SIZE(g_pageStage7Idx), .group=2,.unused={0,0},.params=g_pageStage7Idx },
    { .name="Stage 8",     .numParams=(uint8_t)ARRAY_SIZE(g_pageStage8Idx), .group=2,.unused={0,0},.params=g_pageStage8Idx },
    { .name="Global",      .numParams=(uint8_t)ARRAY_SIZE(g_pageGlobalIdx), .group=0,.unused={0,0},.params=g_pageGlobalIdx },
};
static const _NT_parameterPages g_pages = {
    .numPages = ARRAY_SIZE(g_pages_arr),
    .pages    = g_pages_arr,
};

// ----------------------------
// Runtime helpers
// ----------------------------
static inline float* busPtr(float *busFrames, int busVal, int numFrames) {
    if (!busFrames || busVal <= 0) return nullptr;
    int idx = busVal - 1;
    if (idx < 0 || idx >= kMaxBus) return nullptr;
    return busFrames + idx * numFrames;
}

static float cvLevelToVolts(_MiniMARFA *self, int stage) {
    int raw = self->v[kParamCV1 + stage];
    float u = clampf((float)raw / 1000.0f, 0.0f, 1.0f);
    switch (self->v[kParamVoltageRange]) {
        case 1: return u * 5.0f;           // 0-5V
        case 2: return u * 10.0f - 5.0f;   // +/-5V, NT convenience only
        default:return u * 10.0f;          // 0-10V
    }
}

// Quantize according to voltage range:
//   0-5V  → 1V/oct 12TET semitone grid (1/12 V steps) for pitch sequencing.
//   0-10V → 0.1V steps (MARF original behaviour), suitable for modulation.
//   +/-5V → 0.1V steps, suitable for bipolar modulation.
static float quantizeForRange(_MiniMARFA *self, float v) {
    if (self->v[kParamVoltageRange] == 1)        // 0-5V: pitch
        return roundf(v * 12.0f) / 12.0f;
    return roundf(v * 10.0f) * 0.1f;            // modulation ranges
}

static float entryVoltage(_MiniMARFA *self, int stage) {
    float v = cvLevelToVolts(self, stage);
    if (self->v[kParamQuantCont] == 0)           // Quant
        v = quantizeForRange(self, v);
    return v;
}

static float liveVoltage(_MiniMARFA *self, _MiniMARFA_DTC *d, int stage) {
    if (self->v[kParamQuantCont] == 0) // Quant: sampled at stage entry
        return d->sampledVoltage;
    return cvLevelToVolts(self, stage);
}

static float timeLevelToSeconds(_MiniMARFA *self, int stage) {
    int raw = self->v[kParamTime1 + stage];
    float u = clampf((float)raw / 1000.0f, 0.0f, 1.0f);
    float lo = 0.02f, hi = 0.3f;
    switch (self->v[kParamTimeRange]) {
        case 0: lo = 0.002f; hi = 0.03f; break;
        case 1: lo = 0.02f;  hi = 0.3f;  break;
        case 2: lo = 0.2f;   hi = 3.0f;  break;
        case 3: lo = 2.0f;   hi = 30.0f; break;
    }
    return lo + u * (hi - lo);
}

static int addressedStageFromCV(float v) {
    float vv = clampf(v, 0.0f, 9.999f);
    int s = (int)floorf((vv / 10.0f) * (float)kStages);
    if (s < 0) s = 0;
    if (s >= kStages) s = kStages - 1;
    return s;
}

static void resolveCycleForStage(_MiniMARFA_DTC *d, int stage) {
    bool anyFirst = false;
    bool anyLast = false;
    for (int i=0;i<kStages;i++) {
        if (d->flags[i].first) anyFirst = true;
        if (d->flags[i].last)  anyLast  = true;
    }

    // Manual correction: no FIRST/LAST free-runs the whole range.
    if (!anyFirst && !anyLast) {
        d->cycleFirst = 0;
        d->cycleLast  = kStages - 1;
        return;
    }

    // Find the absolute lowest FIRST stage at or before the current stage,
    // matching strict MARF behaviour where FIRST/LAST are absolute markers.
    int f = 0;
    for (int i=0;i<kStages;i++) {
        if (d->flags[i].first) { f = i; break; }
    }

    int l = kStages - 1;
    for (int i=stage;i<kStages;i++) {
        if (d->flags[i].last) { l = i; break; }
    }

    if (l < f) l = kStages - 1;
    d->cycleFirst = f;
    d->cycleLast  = l;
}

static int nextStage(_MiniMARFA_DTC *d) {
    if (d->currentStage >= d->cycleLast)
        return d->cycleFirst;
    return d->currentStage + 1;
}

static void fireEntryPulses(_MiniMARFA *self, _MiniMARFA_DTC *d) {
    // Pulse Length parameter is in milliseconds; convert to samples.
    int len = (int)roundf((float)self->v[kParamPulseLength] * 0.001f * d->sampleRate);
    if (len < 1) len = 1;

    if (d->flags[d->currentStage].pulse1)
        d->pulse1Remaining = len;
    if (d->flags[d->currentStage].pulse2)
        d->pulse2Remaining = len;

    // MARF all-pulses output fires at the start of every stage.
    d->allPulseRemaining = len;
}

static void enterStage(_MiniMARFA *self, _MiniMARFA_DTC *d, int stage) {
    if (stage < 0) stage = 0;
    if (stage >= kStages) stage = kStages - 1;

    d->previousStage = d->currentStage;
    d->currentStage = stage;
    resolveCycleForStage(d, stage);

    d->phase = 0.0f;
    d->stageDurationSeconds = timeLevelToSeconds(self, stage);
    d->fromVoltage = d->output;
    d->sampledVoltage = entryVoltage(self, stage);
    d->targetVoltage = d->sampledVoltage;

    fireEntryPulses(self, d);

    // Command flags are evaluated at stage entry.
    d->held = false;
    if (d->flags[stage].stop) {
        d->running = false;
        d->held = true;
    } else if (d->flags[stage].sust && d->startGate) {
        d->held = true;
    } else if (d->flags[stage].enable && !d->startGate) {
        d->held = true;
    }
}

static void startAFG(_MiniMARFA *self, _MiniMARFA_DTC *d) {
    if (!d->running) {
        d->running = true;
        // If completely stopped at startup, begin from current stage.
    }

    if (d->held) {
        int s = d->currentStage;
        if (d->flags[s].stop) {
            d->held = false;
            d->running = true;
            enterStage(self, d, nextStage(d));
        } else if (d->flags[s].sust && !d->startGate) {
            d->held = false;
        } else if (d->flags[s].enable && d->startGate) {
            d->held = false;
        }
    }
}

static void resetAFG(_MiniMARFA *self, _MiniMARFA_DTC *d) {
    enterStage(self, d, d->cycleFirst);
}

static void stopAFG(_MiniMARFA_DTC *d) {
    d->running = false;
    d->held = true;
}

static void syncParamsFromSelectedStage(_MiniMARFA *self) {
    // The NT API does not provide a way to push parameter values back to the
    // UI from inside parameterChanged (NT_setParameterValue does not exist).
    // Instead, the flag params (P1/P2/STOP/SUST/ENABLE/FIRST/LAST) in v[]
    // always reflect the most recently *edited* stage, not necessarily the
    // currently *selected* stage.  The draw() function shows the active stage's
    // flags directly from the DTC so the display is always accurate.
    (void)self;
}

static void syncSelectedStageFromParams(_MiniMARFA *self, int p) {
    // Called when any flag param changes.  Determine which stage changed
    // from the param index, then update that stage's flags in the DTC.
    auto *d = self->dtc;
    if (!d) return;
    int s = (p - kParamS1Pulse1) / kFlagsPerStage;
    if (s < 0 || s >= kStages) return;
    int base = kParamS1Pulse1 + s * kFlagsPerStage;

    d->flags[s].pulse1 = (self->v[base + kFlagPulse1] != 0);
    d->flags[s].pulse2 = (self->v[base + kFlagPulse2] != 0);
    d->flags[s].stop   = (self->v[base + kFlagStop]   != 0);
    d->flags[s].sust   = (self->v[base + kFlagSust]   != 0);
    d->flags[s].enable = (self->v[base + kFlagEnable] != 0);
    d->flags[s].first  = (self->v[base + kFlagFirst]  != 0);
    d->flags[s].last   = (self->v[base + kFlagLast]   != 0);
    resolveCycleForStage(d, d->currentStage);
}

static void updateDisplayData(_MiniMARFA *self, _MiniMARFA_DTC *d) {
    // Build a compact time-proportional contour.  The faders already show exact
    // values; this display is for structural/process feedback.
    float times[kStages];
    float total = 0.0f;
    for (int i=0;i<kStages;i++) {
        times[i] = timeLevelToSeconds(self, i);
        total += times[i];
    }
    if (total <= 0.0001f) total = 1.0f;

    int xStage[kStages+1];
    xStage[0] = 0;
    float acc = 0.0f;
    for (int i=0;i<kStages;i++) {
        acc += times[i] / total;
        int x = (int)roundf(acc * (float)(kDisplayPoints-1));
        if (x <= xStage[i]) x = xStage[i] + 1;
        if (x > kDisplayPoints-1) x = kDisplayPoints-1;
        xStage[i+1] = x;
        d->dispStageX[i] = xStage[i];
    }

    for (int i=0;i<kStages;i++) {
        float a = cvLevelToVolts(self, i);
        // Last stage holds flat rather than interpolating back toward stage 0.
        float b = (i < kStages-1) ? cvLevelToVolts(self, i+1) : a;
        int x0 = xStage[i];
        int x1 = xStage[i+1];
        if (x1 <= x0) x1 = x0 + 1;
        for (int x=x0; x<=x1 && x<kDisplayPoints; x++) {
            float t = (float)(x - x0) / (float)(x1 - x0);
            float v = (self->v[kParamSlopedStepped] == 1) ? a : lerpf(a, b, t);
            d->dispY[x] = v;
        }
    }
}

// ----------------------------
// parameterChanged()
// ----------------------------
static void parameterChanged(_NT_algorithm *base, int p) {
    auto *self = (_MiniMARFA*)base;
    auto *d = self->dtc;
    if (!d) return;

    if (p == kParamManualStart)  { startAFG(self, d);  return; }
    if (p == kParamManualStop)   { stopAFG(d);         return; }
    if (p == kParamManualReset)  { resetAFG(self, d);  return; }
    if (p == kParamManualStrobe) {
        // Advance one stage (mirrors the ADVANCE button on the real MARF).
        d->running = true;
        d->held = false;
        enterStage(self, d, nextStage(d));
        return;
    }

    if (p >= kParamS1Pulse1 && p <= (int)kParamS8Last) {
        syncSelectedStageFromParams(self, p);
        return;
    }

    if (p >= kParamCV1 && p <= kParamCV8) {
        if (self->v[kParamQuantCont] == 1 && d->currentStage == (p - kParamCV1))
            d->targetVoltage = cvLevelToVolts(self, d->currentStage);
    }
}

// ----------------------------
// step()
// ----------------------------
static void step(_NT_algorithm *base, float *busFrames, int numFramesBy4) {
    auto *self = (_MiniMARFA*)base;
    auto *d = self->dtc;
    if (!d) return;

    const int numFrames = numFramesBy4 * 4;

    float *vout  = busPtr(busFrames, self->v[kParamVOut],   numFrames);
    float *tout  = busPtr(busFrames, self->v[kParamTOut],   numFrames);
    float *rout  = busPtr(busFrames, self->v[kParamROut],   numFrames);
    float *p1out = busPtr(busFrames, self->v[kParamP1Out],  numFrames);
    float *p2out = busPtr(busFrames, self->v[kParamP2Out],  numFrames);
    float *apout = busPtr(busFrames, self->v[kParamAPOut],  numFrames);

    float *startIn  = busPtr(busFrames, self->v[kParamStartIn],  numFrames);
    float *stopIn   = busPtr(busFrames, self->v[kParamStopIn],   numFrames);
    float *resetIn  = busPtr(busFrames, self->v[kParamResetIn],  numFrames);
    float *strobeIn = busPtr(busFrames, self->v[kParamStrobeIn], numFrames);
    float *sextIn   = busPtr(busFrames, self->v[kParamSExtIn],   numFrames);

    int voutMode = self->v[kParamVOutMode];
    int toutMode = self->v[kParamTOutMode];
    int routMode = self->v[kParamROutMode];
    int p1Mode   = self->v[kParamP1OutMode];
    int p2Mode   = self->v[kParamP2OutMode];
    int apMode   = self->v[kParamAPOutMode];

    for (int n=0;n<numFrames;n++) {
        // Inputs
        if (startIn) {
            bool edge = risingEdge(startIn[n], d->startGate);
            if (edge) startAFG(self, d);
        } else {
            d->startGate = false;
        }
        if (stopIn && risingEdge(stopIn[n], d->stopHigh))
            stopAFG(d);
        if (resetIn && risingEdge(resetIn[n], d->resetHigh))
            resetAFG(self, d);

        if (strobeIn && risingEdge(strobeIn[n], d->strobeHigh)) {
            if (sextIn && self->v[kParamStageAddress] == 1) {
                int target = addressedStageFromCV(sextIn[n]);
                d->running = true;
                d->held = false;
                enterStage(self, d, target);
            }
        }

        if (sextIn && self->v[kParamStageAddress] == 2) {
            int target = addressedStageFromCV(sextIn[n]);
            if (target != d->currentStage) {
                d->running = true;
                d->held = false;
                enterStage(self, d, target);
            }
        }

        // Holds can release as START gate state changes.
        if (d->held) {
            int s = d->currentStage;
            if (d->flags[s].sust && !d->startGate)
                d->held = false;
            if (d->flags[s].enable && d->startGate)
                d->held = false;
        }

        // Output value
        float target = liveVoltage(self, d, d->currentStage);
        float y;
        if (self->v[kParamSlopedStepped] == 1) { // Stepped
            y = target;
        } else { // Sloped
            float t = clampf(d->phase, 0.0f, 1.0f);
            y = lerpf(d->fromVoltage, target, t);
        }
        d->output = y;

        // CV outputs
        if (vout) {
            if (voutMode == 0) vout[n] += y;
            else               vout[n]  = y;
        }
        if (tout) {
            float tv = clampf((float)self->v[kParamTime1 + d->currentStage] / 1000.0f, 0.0f, 1.0f) * 10.0f;
            if (toutMode == 0) tout[n] += tv;
            else               tout[n]  = tv;
        }
        if (rout) {
            float rv = 10.0f * (1.0f - clampf(d->phase, 0.0f, 1.0f));
            if (routMode == 0) rout[n] += rv;
            else               rout[n]  = rv;
        }

        // Pulse outputs
        float p1 = (d->pulse1Remaining > 0) ? kPulseVolts : 0.0f;
        float p2 = (d->pulse2Remaining > 0) ? kPulseVolts : 0.0f;
        float ap = (d->allPulseRemaining > 0) ? kPulseVolts : 0.0f;
        if (p1out) { if (p1Mode == 0) p1out[n] += p1; else p1out[n] = p1; }
        if (p2out) { if (p2Mode == 0) p2out[n] += p2; else p2out[n] = p2; }
        if (apout) { if (apMode == 0) apout[n] += ap; else apout[n] = ap; }

        if (d->pulse1Remaining > 0) --d->pulse1Remaining;
        if (d->pulse2Remaining > 0) --d->pulse2Remaining;
        if (d->allPulseRemaining > 0) --d->allPulseRemaining;

        // Internal timing
        if (d->running && !d->held) {
            float inc = 1.0f / (d->stageDurationSeconds * d->sampleRate);
            d->phase += inc;
            if (d->phase >= 1.0f)
                enterStage(self, d, nextStage(d));
        }
    }

    updateDisplayData(self, d);
}

// ----------------------------
// draw()
// ----------------------------
static bool draw(_NT_algorithm *base) {
    auto *self = (_MiniMARFA*)base;
    auto *d = self->dtc;
    if (!d) return false;

    const int W = 256;
    const int graphTop = 10;
    const int graphBottom = 30;
    const int rowsTop = 34;
    const int rowH = 5;

    auto mapX = [&](int idx)->int {
        return (int)roundf(((float)idx / (float)(kDisplayPoints-1)) * (float)(W-1));
    };
    auto mapY = [&](float v)->int {
        float vv = clampf(v, 0.0f, 10.0f);
        return graphTop + (int)roundf((1.0f - vv / 10.0f) * (float)(graphBottom - graphTop));
    };

    // Contour
    for (int i=1;i<kDisplayPoints;i++) {
        int x0 = mapX(i-1);
        int x1 = mapX(i);
        int y0 = mapY(d->dispY[i-1]);
        int y1 = mapY(d->dispY[i]);
        NT_drawShapeI(kNT_line, x0, y0, x1, y1, 15);
    }

    // Current stage cursor
    int curStart = d->dispStageX[d->currentStage];
    int curEnd   = (d->currentStage < kStages-1) ? d->dispStageX[d->currentStage+1] : kDisplayPoints-1;
    int cx = mapX(curStart + (int)roundf((float)(curEnd - curStart) * clampf(d->phase, 0.0f, 1.0f)));
    NT_drawShapeI(kNT_line, cx, graphTop, cx, graphBottom, 15);

    // Stage boundary ticks
    for (int s=0;s<kStages;s++) {
        int x = mapX(d->dispStageX[s]);
        NT_drawShapeI(kNT_line, x, graphBottom+1, x, graphBottom+4, 7);
    }

    // Programming rows: P1, P2, ST, SU, EN, F, L.
    // Keep labels off-screen-minimal; rows align with stage columns.
    for (int s=0;s<kStages;s++) {
        int x = mapX(d->dispStageX[s]) + 2;
        if (x > W-3) x = W-3;
        bool vals[7] = {
            d->flags[s].pulse1,
            d->flags[s].pulse2,
            d->flags[s].stop,
            d->flags[s].sust,
            d->flags[s].enable,
            d->flags[s].first,
            d->flags[s].last,
        };
        for (int r=0;r<7;r++) {
            if (vals[r]) {
                int y = rowsTop + r * rowH;
                NT_drawShapeI(kNT_rectangle, x, y, x+2, y+2, 15);
            }
        }
    }

    return false;
}

// ----------------------------
// NT registration
// ----------------------------
static void calculateRequirements(_NT_algorithmRequirements &req, const int32_t* specifications) {
    (void)specifications;
    req.numParameters = kNumParams;
    req.sram = sizeof(_MiniMARFA);
    req.dram = 0;
    req.dtc  = sizeof(_MiniMARFA_DTC);
    req.itc  = 0;
}

static _NT_algorithm* construct(const _NT_algorithmMemoryPtrs& ptrs,
                                const _NT_algorithmRequirements& req,
                                const int32_t* specifications) {
    (void)req; (void)specifications;
    auto *d    = new(ptrs.dtc)  _MiniMARFA_DTC();
    auto *self = new(ptrs.sram) _MiniMARFA(d);

    self->parameters     = g_parameters;
    self->parameterPages = &g_pages;

    float sr = (float)NT_globals.sampleRate;
    if (sr < 8000.0f || sr > 192000.0f) sr = kDefaultSampleRate;
    d->sampleRate = sr;

    // Seed DTC flags from parameter values — essential for preset restore,
    // since the NT sets v[] before calling construct.
    for (int s = 0; s < kStages; s++) {
        int base = kParamS1Pulse1 + s * kFlagsPerStage;
        d->flags[s].pulse1 = (self->v[base + kFlagPulse1] != 0);
        d->flags[s].pulse2 = (self->v[base + kFlagPulse2] != 0);
        d->flags[s].stop   = (self->v[base + kFlagStop]   != 0);
        d->flags[s].sust   = (self->v[base + kFlagSust]   != 0);
        d->flags[s].enable = (self->v[base + kFlagEnable] != 0);
        d->flags[s].first  = (self->v[base + kFlagFirst]  != 0);
        d->flags[s].last   = (self->v[base + kFlagLast]   != 0);
    }
    resolveCycleForStage(d, 0);

    // Start at stage 1, held. A START pulse begins traversal.
    d->currentStage = 0;
    d->cycleFirst = 0;
    d->cycleLast = kStages - 1;
    d->sampledVoltage = entryVoltage(self, 0);
    d->targetVoltage = d->sampledVoltage;
    d->output = d->targetVoltage;
    d->fromVoltage = d->output;
    d->stageDurationSeconds = timeLevelToSeconds(self, 0);

    updateDisplayData(self, d);
    return self;
}

static uint32_t hasCustomUi(_NT_algorithm *) { return 0; }
static void customUi(_NT_algorithm *, const _NT_uiData &) {}
static void calculateStaticRequirements(_NT_staticRequirements& req) { req.dram = 0; }
static void initialise(_NT_staticMemoryPtrs& ptrs, const _NT_staticRequirements& req) {
    (void)ptrs; (void)req;
}

static _NT_factory g_factory = {
    .guid        = NT_MULTICHAR('M','M','F','A'), // MiniMARFA internal 4-char ID
    // Gallery GUID/UUID: 2D7699AB-3C80-4D26-9C65-24EED393C4A8
    .name        = "MiniMARFA",
    .description = "8-stage MARF-style arbitrary function programmer",
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
    .tags              = kNT_tagInstrument,
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
