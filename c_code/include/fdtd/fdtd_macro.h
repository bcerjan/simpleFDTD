#ifndef FDTD_MACRO
#define FDTD_MACRO

#include "fdtd_grid.h"

// Main Region Fields:
#define HzG(G) G->hz
#define DahzG(G) G->dahz
#define DbhzG(G) G->dbhz

#define ExG(G) G->ex
#define CaexG(G) G->caex
#define CbexG(G) G->cbex

#define EyG(G) G->ey
#define CaeyG(G) G->caey
#define CbeyG(G) G->cbey

// PML H-field
#define HzyG(G) G->hzy
#define DahzyG(G) G->dahzy
#define DbhzyG(G) G->dbhzy

// Constants:
#define xSizeG(G) G->xSize
#define ySizeG(G) G->ySize
#define maximumIterationG(G) G->maximumIteration
#define xSourceG(G) G->xSource
#define ySourceG(G) G->ySource
#define DxG(G) G->dx
#define DtG(G) G->dt
#define CourantSG(G) G->courantS
#define PiG(G) G->pi
#define SpeedOfLightG(G) G->speedOfLight

// Drude Metals (Polarization Current)
#define JxG(G) G->jx
#define JyG(G) G->jy
#define CjjG(G) G->cjj
#define CjeG(G) G->cje
#define ExOldG(G) G->exOld
#define EyOldG(G) G->eyOld

// Values for tracking DFT:
#define ReflXPosG(G) G->reflXPos
#define TranXPosG(G) G->tranXPos
#define ReflDFTG(G) G->reflDFT
#define TranDFTG(G) G->tranDFT
#define KListG(G) G->kList
#define WavelengthListG(G) G->wavelengthList
#define ReEyDFTG(G) G->reEyDFT
#define ImEyDFTG(G) G->imEyDFT
#define ReHzDFTG(G) G->reHzDFT
#define ImHzDFTG(G) G->imHzDFT

// For plotting:
#define E2FieldG(G) G->e2Field

// Source values and PML Regions
#define SourceValG(G) G->sourceValue
#define RegionDataG(G) G->regionData

#define ObjectLocsG(G) G->object_locs

// Variables for tracking loop state in emscripten
#define TimeStepG(G) G->timeStep
#define IntervalG(G) G->interval

// Now we assume that our grid will be called 'g'

#define hz HzG(g)
#define dahz DahzG(g)
#define dbhz DbhzG(g)

#define ex ExG(g)
#define caex CaexG(g)
#define cbex CbexG(g)

#define ey EyG(g)
#define caey CaeyG(g)
#define cbey CbeyG(g)

// PML H-field
#define hzy HzyG(g)
#define dahzy DahzyG(g)
#define dbhzy DbhzyG(g)

// Constants:
#define xSize xSizeG(g)
#define ySize ySizeG(g)
#define maximumIteration maximumIterationG(g)
#define xSource xSourceG(g)
#define ySource ySourceG(g)
#define dx DxG(g)
#define dt DtG(g)
#define courantS CourantSG(g)
#define pi PiG(g)
#define speedOfLight SpeedOfLightG(g)

// Drude Metals (Polarization Current)
#define jx JxG(g)
#define jy JyG(g)
#define cjj CjjG(g)
#define cje CjeG(g)
#define exOld ExOldG(g)
#define eyOld EyOldG(g)

// Values for tracking DFT:
#define reflXPos ReflXPosG(g)
#define tranXPos TranXPosG(g)
#define reflDFT ReflDFTG(g)
#define tranDFT TranDFTG(g)
#define kList KListG(g)
#define wavelengthList WavelengthListG(g)

#define reEyDFT ReEyDFTG(g)
#define imEyDFT ImEyDFTG(g)
#define reHzDFT ReHzDFTG(g)
#define imHzDFT ImHzDFTG(g)

// For plotting:
#define e2Field E2FieldG(g)

// Source values and PML Regions
#define sourceVal SourceValG(g)
#define regionData RegionDataG(g)

#define object_locs ObjectLocsG(g)

// For emscripten animation:
#define timeStep TimeStepG(g)
#define interval IntervalG(g)

#endif
