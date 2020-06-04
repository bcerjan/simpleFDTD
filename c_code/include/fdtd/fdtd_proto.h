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

#ifndef FDTD_PROTOTYPES
#define FDTD_PROTOTYPES

#include <complex.h>
#include "fdtd_grid.h"

void InitializeFdtd(struct Grid *g, int metalChoice, int objectChoice,
   double objectXSize, double objectYSize,
   double environmentIndex, double objectIndex);
void freeGrid(struct Grid *g);
void freeDoublePtr(double **ptr, int imax);
void freeComplexDoublePtr(complex double **ptr, int imax);
void freeComplexTriplePtr(complex double **ptr, int imax, int jmax);

double **AllocateMemory(int imax, int jmax, double initialValue);
complex double **AllocateComplexMemory(int imax, int jmax, complex double initialValue);
double ***AllocateMemory3D(int imax, int jmax, int kmax, double *initArray);
complex double ***AllocateComplexMemory3D(int imax, int jmax, int kmax, complex double *initArray);
double *AllocateMemory1D(int size, double initialValue);
double *AllocateGridMemory();
void EvaluateFdtd(struct Grid *g, double minimumValue, double maximumValue);

#endif
