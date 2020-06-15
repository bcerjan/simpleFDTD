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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*#include "array_proto.h"
#include "fdtd_proto.h"
#include "fdtd_macro.h"*/
#include "material_data.h"


int main()
{
/*
  struct Grid *g;
  g = AllocateGridMemory();

  InitializeFdtd(g, 0, 0, 100.0, 1.0, 1.0);

  findMatEdge(g);
*/

  printf("params[0] ap: %f + i%f\n", creal(materialData[0].params[0].ap), cimag(materialData[0].params[0].ap));
  printf("params[1] ap: %f + i%f\n", creal(materialData[0].params[1].ap), cimag(materialData[1].params[0].ap));
  printf("conductivity: %f\n", materialData[0].conductivity);
  printf("conversion: %f\n", materialData[0].permeability);
  //freeGrid(g);

  return(0);
}
