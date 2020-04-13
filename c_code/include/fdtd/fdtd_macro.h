/**
    Copyright (c) 2020 Ben Cerjan

    This file is part of simpleFDTD.

    simpleFDTD is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    simpleFDTD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with simpleFDTD.  If not, see <https://www.gnu.org/licenses/>.
**/

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
#define EnvIndexG(G) G->envIndex

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

#define ReEyReflDFTG(G) G->reEyReflDFT
#define ImEyReflDFTG(G) G->imEyReflDFT
#define ReEyTranDFTG(G) G->reEyTranDFT
#define ImEyTranDFTG(G) G->imEyTranDFT

#define ReHzReflDFTG(G) G->reHzReflDFT
#define ImHzReflDFTG(G) G->imHzReflDFT
#define ReHzTranDFTG(G) G->reHzTranDFT
#define ImHzTranDFTG(G) G->imHzTranDFT

// For plotting:
#define E2FieldG(G) G->e2Field
#define EdgeMatG(G) G->edgeMat

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
#define envIndex EnvIndexG(g)

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

#define reEyReflDFT ReEyReflDFTG(g)
#define imEyReflDFT ImEyReflDFTG(g)
#define reEyTranDFT ReEyTranDFTG(g)
#define imEyTranDFT ImEyTranDFTG(g)

#define reHzReflDFT ReHzReflDFTG(g)
#define imHzReflDFT ImHzReflDFTG(g)
#define reHzTranDFT ReHzTranDFTG(g)
#define imHzTranDFT ImHzTranDFTG(g)

// For plotting:
#define e2Field E2FieldG(g)
#define edgeMat EdgeMatG(g)

// Source values and PML Regions
#define sourceVal SourceValG(g)
#define regionData RegionDataG(g)

#define object_locs ObjectLocsG(g)

// For emscripten animation:
#define timeStep TimeStepG(g)
#define interval IntervalG(g)

#endif
