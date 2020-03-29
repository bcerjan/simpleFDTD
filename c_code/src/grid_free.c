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

#include "fdtd_macro.h"
#include <stdlib.h>


void freeDoublePtr(double **ptr, int imax) {
  int i;
  for (i = 0; i < imax; i++) {
    free(ptr[i]);
  } /* iForLoop */

  free(ptr);

  return;
}

void freeGrid(struct Grid *g) {
  // First, free all arrays:
  freeDoublePtr(ex, xSize);
  freeDoublePtr(ey, xSize + 1);
  freeDoublePtr(hz, xSize);
  freeDoublePtr(caex, xSize);
  freeDoublePtr(cbex, xSize);
  freeDoublePtr(caey, xSize);
  freeDoublePtr(cbey, xSize);
  freeDoublePtr(dahz, xSize);
  freeDoublePtr(dbhz, xSize);

  freeDoublePtr(jx, xSize);
  freeDoublePtr(jy, xSize);
  freeDoublePtr(cjj, xSize);
  freeDoublePtr(cje, xSize);
  freeDoublePtr(exOld, xSize);
  freeDoublePtr(eyOld, xSize + 1);
  freeDoublePtr(e2Field, xSize);
  freeDoublePtr(edgeMat, xSize);
  freeDoublePtr(object_locs, xSize);

  freeDoublePtr(reEyReflDFT, NUMBERDFTFREQS);
  freeDoublePtr(imEyReflDFT, NUMBERDFTFREQS);
  freeDoublePtr(reEyTranDFT, NUMBERDFTFREQS);
  freeDoublePtr(imEyTranDFT, NUMBERDFTFREQS);
  freeDoublePtr(reHzReflDFT, NUMBERDFTFREQS);
  freeDoublePtr(imHzReflDFT, NUMBERDFTFREQS);
  freeDoublePtr(reHzTranDFT, NUMBERDFTFREQS);
  freeDoublePtr(imHzTranDFT, NUMBERDFTFREQS);

  // Now, all individual pointers:
  free(hzy);
  free(dahzy);
  free(dbhzy);
  free(reflDFT);
  free(tranDFT);

  // And finally the grid itself:
  free(g);

  return;
}
