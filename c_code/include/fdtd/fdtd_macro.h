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
#define ExG(G) G->ex
#define EyG(G) G->ey
#define ExOldG(G) G->exOld
#define EyOldG(G) G->eyOld
#define HzOldG(G) G->hzOld


#define QxG(G) G->qx
#define QyG(G) G->qy
#define QSumCG(G) G->qSumC

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

// Update Equation Constants:
#define HeConstG(G) G->heConst
#define EhConstG(G) G->ehConst
#define EqConstG(G) G->eqConst
#define ABConstG(G) G->ABConst
#define QConst1G(G) G->qConst1
#define QConst2G(G) G->qConst2
#define QxSumG(G) G->qxSum
#define QySumG(G) G->qySum
#define IConst1G(G) G->iConst1
#define IConst2G(G) G->iConst2
#define Number_polesG(G) G->number_poles

// Tridiagonal values:
#define AexG(G) G->aex
#define BexG(G) G->bex
#define CexG(G) G->cex
#define AhzG(G) G->ahz
#define BhzG(G) G->bhz
#define ChzG(G) G->chz

// Absorbing boundary factor:
#define AbsConstG(G) G->absConst


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
#define RefractiveIndexIndexG(G) G->refractiveIndexIndex

// Variables for tracking loop state in emscripten
#define TimeStepG(G) G->timeStep
#define IntervalG(G) G->interval

// Now we assume that our grid will be called 'g'

#define hz HzG(g)
#define ex ExG(g)
#define ey EyG(g)

#define exOld ExOldG(g)
#define eyOld EyOldG(g)
#define hzOld HzOldG(g)

#define qx QxG(g)
#define qy QyG(g)
#define qSumC QSumCG(g)

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

// Update Equation Constants:
#define heConst HeConstG(g)
#define ehConst EhConstG(g)
#define eqConst EqConstG(g)
#define ABConst ABConstG(g)
#define qConst1 QConst1G(g)
#define qConst2 QConst2G(g)
#define qxSum QxSumG(g)
#define qySum QySumG(g)
#define iConst1 IConst1G(g)
#define iConst2 IConst2G(g)
#define number_poles Number_polesG(g)
// Tridiagonal values:
#define aex AexG(g)
#define bex BexG(g)
#define cex CexG(g)
#define ahz AhzG(g)
#define bhz BhzG(g)
#define chz ChzG(g)

#define absConst AbsConstG(g)

#define PMLkx PMLkxG(g)
#define PMLky PMLkyG(g)

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
#define refractiveIndexIndex RefractiveIndexIndexG(g)

// For emscripten animation:
#define timeStep TimeStepG(g)
#define interval IntervalG(g)

#endif
