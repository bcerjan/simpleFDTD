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

#ifndef STRUCT_FUNCS
#define STRUCT_FUNCS

#include "fdtd_macro.h"

void structInit(int xCenter, int yCenter);
void addDisk(struct Grid *g, double xRadius, double yRadius);
void addRect(struct Grid *g, double width, double length);
void addTriangle(struct Grid *g, double xLength, double yLength);

#endif
