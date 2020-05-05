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

#define EyG(G) G->ey

// PML H-field
#define HGrad1G(G) G->hGrad1
#define EGrad1G(G) G->eGrad1
#define HGrad2G(G) G->hGrad2
#define EGrad2G(G) G->eGrad2
#define HGrad3G(G) G->hGrad3
#define EGrad3G(G) G->eGrad3

#define PmlSxG(G) G->pmlSx
#define PmlSyG(G) G->pmlSy
#define PmlTzG(G) G->PmlTz
#define PmlSxOldG(G) G->pmlSxOld
#define PmlSyOldG(G) G->pmlSyOld
#define PmlTzOldG(G) G->PmlTzOld

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

// Critical-Point Material Data (Polarization)
#define Number_polesG(G) G->number_poles
#define PxG(G) G->px
#define PyG(G) G->py
#define PxOldG(G) G->pxOld
#define PyOldG(G) G->pyOld
#define ExOldG(G) G->exOld
#define EyOldG(G) G->eyOld
#define ExOld2G(G) G->exOld2
#define EyOld2G(G) G->eyOld2
#define HzOldG(G) G->hzOld

#define C1SumXG(G) G->c1SumX
#define C2SumXG(G) G->c2SumX
#define C1SumYG(G) G->c1SumY
#define C2SumYG(G) G->c2SumY
#define C3SumG(G) G->c3Sum
#define C4SumG(G) G->c4Sum
#define C5SumG(G) G->c5Sum
#define C1GridG(G) G->c1Grid
#define C2GridG(G) G->c2Grid
#define C3GridG(G) G->c3Grid
#define C4GridG(G) G->c4Grid
#define C5GridG(G) G->c5Grid

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
#define dahz DahzG(g)
#define dbhz DbhzG(g)

#define ex ExG(g)

#define ey EyG(g)

// PML Coefficients
#define hGrad1 HGrad1G(g)
#define eGrad1 EGrad1G(g)
#define hGrad2 HGrad2G(g)
#define eGrad2 EGrad2G(g)
#define hGrad3 HGrad3G(g)
#define eGrad3 EGrad3G(g)

#define pmlSx PmlSxG(g)
#define pmlSy PmlSyG(g)
#define pmlTz PmlTzG(g)
#define pmlSxOld PmlSxOldG(g)
#define pmlSyOld PmlSyOldG(g)
#define pmlTzOld PmlTzOldG(g)

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

// Critical Point Materials (Polarization)
#define number_poles Number_polesG(g)
#define px PxG(g)
#define py PyG(g)
#define pxOld PxOldG(g)
#define pyOld PyOldG(g)
#define exOld ExOldG(g)
#define eyOld EyOldG(g)
#define exOld2 ExOld2G(g)
#define eyOld2 EyOld2G(g)
#define hzOld HzOldG(g)

#define c1SumX C1SumXG(g)
#define c2SumX C2SumXG(g)
#define c1SumY C1SumYG(g)
#define c2SumY C2SumYG(g)
#define c3Sum C3SumG(g)
#define c4Sum C4SumG(g)
#define c5Sum C5SumG(g)
#define c1Grid C1GridG(g)
#define c2Grid C2GridG(g)
#define c3Grid C3GridG(g)
#define c4Grid C4GridG(g)
#define c5Grid C5GridG(g)

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
