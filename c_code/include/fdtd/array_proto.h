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

#ifndef ARRAYFUNCS
#define ARRAYFUNCS

#include "fdtd_grid.h"

double ArrayMax(double **ptr, int xWidth, int yWidth);
double ArrayMin(double **ptr, int xWidth, int yWidth);
double AbsArrayMax(double **ptr, int xWidth, int yWidth);
double AbsArrayMin(double **ptr, int xWidth, int yWidth);
double AbsVectorMax(double *ptr, int length);

#endif
