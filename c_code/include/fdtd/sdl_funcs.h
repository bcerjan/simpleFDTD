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

#ifndef SDL_IMAGE_FUNC
#define SDL_IMAGE_FUNC

#include "fdtd_macro.h"

void imageInit(struct Grid *g);
void imageFree();
void findMatEdge(struct Grid *g);
void imageShow(struct Grid *g);
void PlotField (struct Grid *g, double maximumValue, double minimumValue);

#endif
