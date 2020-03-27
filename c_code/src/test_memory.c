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
#include "fdtd_proto.h"
#include "fdtd_macro.h"

int main()
{

  int m,n;
  m = 500;
  n = 501;

  struct Grid *g;
  g = AllocateGridMemory();

  InitializeFdtd(g, 0);

/*  double **ptr1,**ptr2,**ptr3, **ptr4, **ptr5, **ptr6;
  ptr1 = AllocateMemory(m, n, 0.0);
  ptr2 = AllocateMemory(m, n, 0.0);
  ptr3 = AllocateMemory(m, n, 0.0);
  ptr4 = AllocateMemory(m, n, 0.0);
  ptr5 = AllocateMemory(m, n, 0.0);
  ptr6 = AllocateMemory(m, n, 0.0);

  ex = AllocateMemory(m, n, 1.0);
  ey = AllocateMemory(m, n, 1.0);
  hz = AllocateMemory(m, n, 1.0);
*/
  printf("ex[0][0]: %f\n", ex[0][0]);
  freeGrid(g);

  return(0);
}
