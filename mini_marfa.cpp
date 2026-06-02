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
#include <cstdio>

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
    bool  pulse1;
    bool  pulse2;
    bool  stop;
    bool  sust;
    bool  enable;
    bool  first;
    bool  last;
    bool  slope;
    bool  vext;
    int   shape;   // -100 to +100
    int   octave;  // 0-3

    StageFlags()
    : pulse1(false), pulse2(false), stop(false), sust(false),
      enable(false), first(false), last(false),
      slope(false), vext(false), shape(0), octave(0) {}
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

    // Set by parameterChanged for manual triggers, consumed in step()
    bool  manualStart;
    bool  manualStop;
    bool  manualReset;
    bool  manualStrobe;

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
        manualStart = manualStop = manualReset = manualStrobe = false;

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
    kParamS8Last   = kParamS1Pulse1 + 8*11 - 1,

    // Global page — starts immediately after the 56 stage flag params
    kParamManualStart  = kParamS8Last + 1,
    kParamManualStop,
    kParamManualReset,
    kParamManualStrobe,
    kParamClearProg,
    kParamQuantCont,
    kParamTimeRange,
    kParamVoltageRange,
    kParamStageAddress,
    kParamPulseLength,
    kParamScale,
    kParamRootNote,
    kParamTimeMult,
    kParamTimeMultIn,
    kParamExtCVIn,
};
static constexpr int kNumParams = kParamExtCVIn + 1;
static constexpr int kFlagsPerStage = 11;
static constexpr int kFlagPulse1 = 0, kFlagPulse2 = 1, kFlagStop  = 2,
                     kFlagSust   = 3, kFlagEnable = 4, kFlagFirst = 5,
                     kFlagLast   = 6, kFlagCurve  = 7, kFlagSlope = 8,
                     kFlagVExt   = 9, kFlagOctave = 10;

// ----------------------------
// String tables
// ----------------------------
static const char* onOffStrings[]          = {"Off", "On", nullptr};
static const char* quantContStrings[]      = {"Quant", "Cont", nullptr};
static const char* timeRangeStrings[]      = {".002-.03s", ".02-.3s", ".2-3s", "2-30s", nullptr};
static const char* voltageRangeStrings[]   = {"0-10V", "0-5V", "0-2V", "+/-5V", nullptr};
static const char* octaveStrings[]         = {"+0", "+1", "+2", "+3", nullptr};
static const char* scaleStrings[]          = {"Off","Major","Minor","Pent Maj","Pent Min", nullptr};
static const char* rootStrings[]           = {"C","C#","D","D#","E","F","F#","G","G#","A","A#","B", nullptr};
static const char* stageAddressStrings[]   = {"Internal", "Strobe Ext", "Cont Ext", nullptr};

// ----------------------------
// Parameter descriptors
// ----------------------------
static const _NT_parameter g_parameters[kNumParams] = {
    // Routing. Outputs use WITH_MODE to match the working Cessna pattern.
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Function CV Out", 1, 13)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Time Out",        0, 0)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Ref Out",         0, 0)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Pulse 1 Out",     0, 0)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Pulse 2 Out",     0, 0)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("All Pulses Out",  0, 0)
    NT_PARAMETER_AUDIO_INPUT("Start Input",     0, 0)
    NT_PARAMETER_AUDIO_INPUT("Stop Input",      0, 0)
    NT_PARAMETER_AUDIO_INPUT("Reset Input",     0, 0)
    NT_PARAMETER_AUDIO_INPUT("Strobe Input",    0, 0)
    NT_PARAMETER_AUDIO_INPUT("Stage Ext Input", 0, 0)

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
    { .name="Last",   .min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings }, \
    { .name="Shape",  .min=-100,.max=100,.def=0,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr }, \
    { .name="Slope",  .min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings }, \
    { .name="V.Ext",  .min=0,.max=1,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=onOffStrings }, \
    { .name="Octave", .min=0,.max=3,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=octaveStrings },
    STAGE_FLAGS STAGE_FLAGS STAGE_FLAGS STAGE_FLAGS
    STAGE_FLAGS STAGE_FLAGS STAGE_FLAGS STAGE_FLAGS
#undef STAGE_FLAGS


    // Global — manual trigger buttons.
    // min=-1/max=1/def=0: encoder can always move either direction, fires every time.
    { .name="Start",      .min=-1,.max=1,.def=0,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Stop",       .min=-1,.max=1,.def=0,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Reset",      .min=-1,.max=1,.def=0,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Strobe",     .min=-1,.max=1,.def=0,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Clear Prog", .min=-1,.max=1,.def=0,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Quant/Cont",    .min=0,.max=1,.def=1,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=quantContStrings },
    { .name="Time Range",    .min=0,.max=3,.def=1,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=timeRangeStrings },
    { .name="Voltage Range", .min=0,.max=3,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=voltageRangeStrings },
    { .name="Stage Address", .min=0,.max=2,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=stageAddressStrings },
    { .name="Pulse ms",      .min=1,.max=100,.def=1,.unit=kNT_unitMs,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    { .name="Scale",         .min=0,.max=4,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=scaleStrings },
    { .name="Root Note",     .min=0,.max=11,.def=0,.unit=kNT_unitEnum,.scaling=kNT_scalingNone,.enumStrings=rootStrings },
    { .name="Time Mult",     .min=1,.max=200,.def=100,.unit=kNT_unitNone,.scaling=kNT_scalingNone,.enumStrings=nullptr },
    NT_PARAMETER_CV_INPUT("Time Mult In", 0, 0)
    NT_PARAMETER_CV_INPUT("Ext CV In",    0, 0)
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
    (uint8_t)kParamResetIn,(uint8_t)kParamStrobeIn,
    (uint8_t)kParamSExtIn,(uint8_t)kParamTimeMultIn,(uint8_t)kParamExtCVIn,
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
    (uint8_t)(kParamS1Pulse1+6),(uint8_t)(kParamS1Pulse1+7),(uint8_t)(kParamS1Pulse1+8),
    (uint8_t)(kParamS1Pulse1+9),(uint8_t)(kParamS1Pulse1+10),
};
static const uint8_t g_pageStage2Idx[] = {
    (uint8_t)(kParamS1Pulse1+11),(uint8_t)(kParamS1Pulse1+12),(uint8_t)(kParamS1Pulse1+13),
    (uint8_t)(kParamS1Pulse1+14),(uint8_t)(kParamS1Pulse1+15),(uint8_t)(kParamS1Pulse1+16),
    (uint8_t)(kParamS1Pulse1+17),(uint8_t)(kParamS1Pulse1+18),(uint8_t)(kParamS1Pulse1+19),
    (uint8_t)(kParamS1Pulse1+20),(uint8_t)(kParamS1Pulse1+21),
};
static const uint8_t g_pageStage3Idx[] = {
    (uint8_t)(kParamS1Pulse1+22),(uint8_t)(kParamS1Pulse1+23),(uint8_t)(kParamS1Pulse1+24),
    (uint8_t)(kParamS1Pulse1+25),(uint8_t)(kParamS1Pulse1+26),(uint8_t)(kParamS1Pulse1+27),
    (uint8_t)(kParamS1Pulse1+28),(uint8_t)(kParamS1Pulse1+29),(uint8_t)(kParamS1Pulse1+30),
    (uint8_t)(kParamS1Pulse1+31),(uint8_t)(kParamS1Pulse1+32),
};
static const uint8_t g_pageStage4Idx[] = {
    (uint8_t)(kParamS1Pulse1+33),(uint8_t)(kParamS1Pulse1+34),(uint8_t)(kParamS1Pulse1+35),
    (uint8_t)(kParamS1Pulse1+36),(uint8_t)(kParamS1Pulse1+37),(uint8_t)(kParamS1Pulse1+38),
    (uint8_t)(kParamS1Pulse1+39),(uint8_t)(kParamS1Pulse1+40),(uint8_t)(kParamS1Pulse1+41),
    (uint8_t)(kParamS1Pulse1+42),(uint8_t)(kParamS1Pulse1+43),
};
static const uint8_t g_pageStage5Idx[] = {
    (uint8_t)(kParamS1Pulse1+44),(uint8_t)(kParamS1Pulse1+45),(uint8_t)(kParamS1Pulse1+46),
    (uint8_t)(kParamS1Pulse1+47),(uint8_t)(kParamS1Pulse1+48),(uint8_t)(kParamS1Pulse1+49),
    (uint8_t)(kParamS1Pulse1+50),(uint8_t)(kParamS1Pulse1+51),(uint8_t)(kParamS1Pulse1+52),
    (uint8_t)(kParamS1Pulse1+53),(uint8_t)(kParamS1Pulse1+54),
};
static const uint8_t g_pageStage6Idx[] = {
    (uint8_t)(kParamS1Pulse1+55),(uint8_t)(kParamS1Pulse1+56),(uint8_t)(kParamS1Pulse1+57),
    (uint8_t)(kParamS1Pulse1+58),(uint8_t)(kParamS1Pulse1+59),(uint8_t)(kParamS1Pulse1+60),
    (uint8_t)(kParamS1Pulse1+61),(uint8_t)(kParamS1Pulse1+62),(uint8_t)(kParamS1Pulse1+63),
    (uint8_t)(kParamS1Pulse1+64),(uint8_t)(kParamS1Pulse1+65),
};
static const uint8_t g_pageStage7Idx[] = {
    (uint8_t)(kParamS1Pulse1+66),(uint8_t)(kParamS1Pulse1+67),(uint8_t)(kParamS1Pulse1+68),
    (uint8_t)(kParamS1Pulse1+69),(uint8_t)(kParamS1Pulse1+70),(uint8_t)(kParamS1Pulse1+71),
    (uint8_t)(kParamS1Pulse1+72),(uint8_t)(kParamS1Pulse1+73),(uint8_t)(kParamS1Pulse1+74),
    (uint8_t)(kParamS1Pulse1+75),(uint8_t)(kParamS1Pulse1+76),
};
static const uint8_t g_pageStage8Idx[] = {
    (uint8_t)(kParamS1Pulse1+77),(uint8_t)(kParamS1Pulse1+78),(uint8_t)(kParamS1Pulse1+79),
    (uint8_t)(kParamS1Pulse1+80),(uint8_t)(kParamS1Pulse1+81),(uint8_t)(kParamS1Pulse1+82),
    (uint8_t)(kParamS1Pulse1+83),(uint8_t)(kParamS1Pulse1+84),(uint8_t)(kParamS1Pulse1+85),
    (uint8_t)(kParamS1Pulse1+86),(uint8_t)(kParamS1Pulse1+87),
};
static const uint8_t g_pageGlobalIdx[] = {
    (uint8_t)kParamManualStart,
    (uint8_t)kParamManualStop,
    (uint8_t)kParamManualReset,
    (uint8_t)kParamManualStrobe,
    (uint8_t)kParamClearProg,
    (uint8_t)kParamQuantCont,(uint8_t)kParamTimeRange,
    (uint8_t)kParamVoltageRange,(uint8_t)kParamStageAddress,(uint8_t)kParamPulseLength,
    (uint8_t)kParamScale,(uint8_t)kParamRootNote,
    (uint8_t)kParamTimeMult,
};

static const _NT_parameterPage g_pages_arr[] = {
    { .name="Global",      .numParams=(uint8_t)ARRAY_SIZE(g_pageGlobalIdx), .group=0,.unused={0,0},.params=g_pageGlobalIdx },
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
        case 2: return u * 2.0f;           // 0-2V: 1V/oct over 2 octaves
        case 3: return u * 10.0f - 5.0f;   // +/-5V
        default:return u * 10.0f;          // 0-10V
    }
}

// Scale intervals in semitones from root
static const int kScaleMajor[]    = {0,2,4,5,7,9,11};
static const int kScaleMinor[]    = {0,2,3,5,7,8,10};
static const int kScalePentMaj[]  = {0,2,4,7,9};
static const int kScalePentMin[]  = {0,3,5,7,10};
static const int kScaleLens[]     = {7,7,5,5};
static const int* kScales[]       = {kScaleMajor,kScaleMinor,kScalePentMaj,kScalePentMin};

static float quantizeToScale(float v, int scaleIdx, int root) {
    // Convert voltage to semitones (1V/oct)
    float semitones = v * 12.0f;
    int octave      = (int)floorf(semitones / 12.0f);
    float semInOct  = semitones - octave * 12.0f;

    // Find nearest scale degree relative to root
    const int* scale = kScales[scaleIdx];
    int len          = kScaleLens[scaleIdx];
    int best         = 0;
    float bestDist   = 999.0f;
    for (int i = 0; i < len; i++) {
        float deg  = fmodf((float)((scale[i] + root) % 12), 12.0f);
        float dist = fabsf(semInOct - deg);
        // Also check wrapping distance
        if (dist > 6.0f) dist = 12.0f - dist;
        if (dist < bestDist) { bestDist = dist; best = (scale[i] + root) % 12; }
    }
    return (float)(octave * 12 + best) / 12.0f;
}

// Quantize according to voltage range and scale settings.
//   0-2V / 0-5V → 1V/oct pitch grid, with optional scale quantization
//   0-10V / +/-5V → 0.1V steps (modulation), scale ignored
static float quantizeForRange(_MiniMARFA *self, float v) {
    int range = self->v[kParamVoltageRange];
    if (range == 1 || range == 2) {  // pitch ranges
        int scale = self->v[kParamScale];
        if (scale > 0)
            return quantizeToScale(v, scale - 1, self->v[kParamRootNote]);
        return roundf(v * 12.0f) / 12.0f;  // chromatic
    }
    return roundf(v * 10.0f) * 0.1f;  // modulation
}

// Returns the base voltage for a stage, applying VEXT and octave offset.
// extCV: current value of the Ext CV Input (pass 0.0f if not available at entry time).
static float stageVoltage(_MiniMARFA *self, _MiniMARFA_DTC *d, int stage, float extCV) {
    bool vext   = d->flags[stage].vext;
    int  octave = d->flags[stage].octave;
    float v = vext ? extCV : cvLevelToVolts(self, stage);
    v += (float)octave;
    return v;
}

static float entryVoltage(_MiniMARFA *self, _MiniMARFA_DTC *d, int stage, float extCV) {
    float v = stageVoltage(self, d, stage, extCV);
    if (self->v[kParamQuantCont] == 0)
        v = quantizeForRange(self, v);
    return v;
}

static float liveVoltage(_MiniMARFA *self, _MiniMARFA_DTC *d, int stage, float extCV) {
    if (self->v[kParamQuantCont] == 0)
        return d->sampledVoltage;
    return stageVoltage(self, d, stage, extCV);
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

static void enterStage(_MiniMARFA *self, _MiniMARFA_DTC *d, int stage, float extCV = 0.0f) {
    if (stage < 0) stage = 0;
    if (stage >= kStages) stage = kStages - 1;

    d->previousStage = d->currentStage;
    d->currentStage = stage;
    resolveCycleForStage(d, stage);

    d->phase = 0.0f;
    d->stageDurationSeconds = timeLevelToSeconds(self, stage);
    d->fromVoltage = d->output;
    d->sampledVoltage = entryVoltage(self, d, stage, extCV);
    d->targetVoltage = d->sampledVoltage;

    fireEntryPulses(self, d);

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
    if (d->held) {
        int s = d->currentStage;
        if (d->flags[s].stop) {
            // Held at a STOP stage — advance to next stage.
            d->held = false;
            d->running = true;
            enterStage(self, d, nextStage(d));
        } else if (d->flags[s].sust && !d->startGate) {
            // Held at SUST with gate low — release hold, continue.
            d->held = false;
        } else if (d->flags[s].enable && d->startGate) {
            // Held at ENABLE with gate high — release hold, continue.
            d->held = false;
        } else {
            // Generic held state (e.g. initial startup) — just start from cycleFirst.
            d->held = false;
            d->running = true;
            enterStage(self, d, d->cycleFirst);
        }
        return;
    }

    if (!d->running) {
        d->running = true;
        enterStage(self, d, d->cycleFirst);
    }
}

static void resetAFG(_MiniMARFA *self, _MiniMARFA_DTC *d) {
    enterStage(self, d, d->cycleFirst);
}

static void stopAFG(_MiniMARFA_DTC *d) {
    d->running = false;
    d->held = true;
}

static void clearProg(_MiniMARFA_DTC *d) {
    for (int s = 0; s < kStages; s++)
        d->flags[s] = {};
    resolveCycleForStage(d, d->currentStage);
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
    d->flags[s].shape  =  self->v[base + kFlagCurve];
    d->flags[s].slope  = (self->v[base + kFlagSlope]  != 0);
    d->flags[s].vext   = (self->v[base + kFlagVExt]   != 0);
    d->flags[s].octave =  self->v[base + kFlagOctave];
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
        float b = (i < kStages-1) ? cvLevelToVolts(self, i+1) : a;
        int x0 = xStage[i];
        int x1 = xStage[i+1];
        if (x1 <= x0) x1 = x0 + 1;
        bool sloped = d->flags[i].slope;
        int shape   = d->flags[i].shape;
        for (int x=x0; x<=x1 && x<kDisplayPoints; x++) {
            float t = (float)(x - x0) / (float)(x1 - x0);
            if (sloped && shape != 0) {
                float exponent = powf(3.0f, (float)shape / 100.0f);
                t = powf(t, exponent);
            }
            float v = sloped ? lerpf(a, b, t) : a;
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

    // Trigger params: set a flag in DTC, consumed safely in step()
    if (p == kParamManualStart)  { d->manualStart  = true; return; }
    if (p == kParamManualStop)   { d->manualStop   = true; return; }
    if (p == kParamManualReset)  { d->manualReset  = true; return; }
    if (p == kParamManualStrobe) { d->manualStrobe = true; return; }
    if (p == kParamClearProg)    { clearProg(d); return; }

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
    float *timeMultIn = busPtr(busFrames, self->v[kParamTimeMultIn], numFrames);
    float *extCVIn    = busPtr(busFrames, self->v[kParamExtCVIn],   numFrames);

    int voutMode = self->v[kParamVOutMode];
    int toutMode = self->v[kParamTOutMode];
    int routMode = self->v[kParamROutMode];
    int p1Mode   = self->v[kParamP1OutMode];
    int p2Mode   = self->v[kParamP2OutMode];
    int apMode   = self->v[kParamAPOutMode];

    for (int n=0;n<numFrames;n++) {
        float extCV = extCVIn ? extCVIn[n] : 0.0f;

        // Manual triggers from parameterChanged (consumed once per block start)
        if (n == 0) {
            if (d->manualStart)  { d->manualStart  = false; startAFG(self, d); }
            if (d->manualStop)   { d->manualStop   = false; stopAFG(d); }
            if (d->manualReset)  { d->manualReset  = false; resetAFG(self, d); }
            if (d->manualStrobe) { d->manualStrobe = false;
                                   d->running = true; d->held = false;
                                   enterStage(self, d, nextStage(d), extCV); }
        }

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
                int tgt = addressedStageFromCV(sextIn[n]);
                d->running = true;
                d->held = false;
                enterStage(self, d, tgt, extCV);
            }
        }

        if (sextIn && self->v[kParamStageAddress] == 2) {
            int tgt = addressedStageFromCV(sextIn[n]);
            if (tgt != d->currentStage) {
                d->running = true;
                d->held = false;
                enterStage(self, d, tgt, extCV);
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
        float target = liveVoltage(self, d, d->currentStage, extCV);
        float y;
        bool sloped = d->flags[d->currentStage].slope;
        if (!sloped) { // Stepped
            y = target;
        } else { // Sloped
            float t = clampf(d->phase, 0.0f, 1.0f);
            // Apply per-stage shape from DTC (-100=log, 0=lin, +100=exp)
            int shape = d->flags[d->currentStage].shape;
            if (shape != 0) {
                float exponent = powf(3.0f, (float)shape / 100.0f);
                t = powf(t, exponent);
            }
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
            // Time multiplier: param is 1-200 representing 0.01x-2.0x
            // CV input adds 0-10V mapped to 0-1.0x additional multiplier
            float mult = (float)self->v[kParamTimeMult] / 100.0f;
            if (timeMultIn) mult += clampf(timeMultIn[n] / 10.0f, 0.0f, 1.0f);
            mult = clampf(mult, 0.01f, 4.0f);  // safety clamp
            float inc = mult / (d->stageDurationSeconds * d->sampleRate);
            d->phase += inc;
            if (d->phase >= 1.0f)
                enterStage(self, d, nextStage(d), extCV);
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

    const int W           = 256;
    const int graphTop    = 10;              // below NT parameter label (~8px)
    const int graphBottom = 42;              // 32px for graph
    const int rowsTop     = graphBottom + 7; // 5px more gap below graph

    auto mapX = [&](int idx)->int {
        return (int)roundf(((float)idx / (float)(kDisplayPoints-1)) * (float)(W-1));
    };
    // Determine display voltage range based on current setting
    float vMin, vMax;
    switch (self->v[kParamVoltageRange]) {
        case 1: vMin = 0.0f; vMax = 5.0f;  break;  // 0-5V
        case 2: vMin = 0.0f; vMax = 2.0f;  break;  // 0-2V
        case 3: vMin =-5.0f; vMax = 5.0f;  break;  // +/-5V
        default:vMin = 0.0f; vMax = 10.0f; break;  // 0-10V
    }

    auto mapY = [&](float v)->int {
        float vv = clampf(v, vMin, vMax);
        return graphTop + (int)roundf((1.0f - (vv - vMin) / (vMax - vMin)) * (float)(graphBottom - graphTop));
    };

    // Contour
    for (int i=1;i<kDisplayPoints;i++) {
        int x0 = mapX(i-1), x1 = mapX(i);
        int y0 = mapY(d->dispY[i-1]), y1 = mapY(d->dispY[i]);
        NT_drawShapeI(kNT_line, x0, y0, x1, y1, 15);
    }

    // Current stage cursor
    int curStart = d->dispStageX[d->currentStage];
    int curEnd   = (d->currentStage < kStages-1) ? d->dispStageX[d->currentStage+1] : kDisplayPoints-1;
    int cx = mapX(curStart + (int)roundf((float)(curEnd - curStart) * clampf(d->phase, 0.0f, 1.0f)));
    NT_drawShapeI(kNT_line, cx, graphTop, cx, graphBottom, 15);

    // Flag display: 2 rows of equal-width columns, independent of stage time.
    // Row 0: P1 P2 ST  |  Row 1: SU EN F L
    const int colW  = W / kStages;   // 32px per stage
    const int rowH2 = 8;
    const int row0y   = rowsTop;
    const int row1y   = rowsTop + rowH2 + 1;

    // Legend: fixed labels on left of each row
    int legendY = row0y - 4;
    NT_drawText(0,   legendY, "Row1: P1 P2 STP", 7, kNT_textLeft,  kNT_textTiny);
    NT_drawText(255, legendY, "Row2: SUS ENA FST LST", 7, kNT_textRight, kNT_textTiny);

    for (int s = 0; s < kStages; s++) {
        int x = s * colW;
        NT_drawShapeI(kNT_line, x, row0y, x, row1y + rowH2, 5);

        int t = colW / 3;
        if (d->flags[s].pulse1) NT_drawShapeI(kNT_rectangle, x+1,      row0y+1, x+t-1,    row0y+rowH2-1, 15);
        if (d->flags[s].pulse2) NT_drawShapeI(kNT_rectangle, x+t+1,    row0y+1, x+2*t-1,  row0y+rowH2-1, 15);
        if (d->flags[s].stop)   NT_drawShapeI(kNT_rectangle, x+2*t+1,  row0y+1, x+colW-2, row0y+rowH2-1, 15);

        int q = colW / 4;
        if (d->flags[s].sust)   NT_drawShapeI(kNT_rectangle, x+1,      row1y+1, x+q-1,    row1y+rowH2-1, 15);
        if (d->flags[s].enable) NT_drawShapeI(kNT_rectangle, x+q+1,    row1y+1, x+2*q-1,  row1y+rowH2-1, 15);
        if (d->flags[s].first)  NT_drawShapeI(kNT_rectangle, x+2*q+1,  row1y+1, x+3*q-1,  row1y+rowH2-1, 15);
        if (d->flags[s].last)   NT_drawShapeI(kNT_rectangle, x+3*q+1,  row1y+1, x+colW-2, row1y+rowH2-1, 15);
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
        d->flags[s].shape  =  self->v[base + kFlagCurve];
        d->flags[s].slope  = (self->v[base + kFlagSlope]  != 0);
        d->flags[s].vext   = (self->v[base + kFlagVExt]   != 0);
        d->flags[s].octave =  self->v[base + kFlagOctave];
    }
    resolveCycleForStage(d, 0);

    // Start halted at cycleFirst. A START pulse begins traversal.
    d->currentStage   = d->cycleFirst;
    d->sampledVoltage = entryVoltage(self, d, d->cycleFirst, 0.0f);
    d->targetVoltage  = d->sampledVoltage;
    d->output         = d->targetVoltage;
    d->fromVoltage    = d->output;
    d->stageDurationSeconds = timeLevelToSeconds(self, d->cycleFirst);

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
