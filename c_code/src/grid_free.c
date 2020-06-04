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
#include <complex.h>


void freeDoublePtr(double **ptr, int imax) {
  int i;
  for (i = 0; i < imax; i++) {
    free(ptr[i]);
  } /* iForLoop */

  free(ptr);

  return;
}

void freeTriplePtr(double ***ptr, int imax, int jmax) {
  int i;
  for (i = 0; i < imax; i++) {
    freeDoublePtr(ptr[i],jmax);
  }
  free(ptr);
  return;
}

void freeComplexDoublePtr(complex double **ptr, int imax) {
  int i;
  for (i = 0; i < imax; i++) {
    free(ptr[i]);
  } /* iForLoop */

  free(ptr);

  return;
}

void freeComplexTriplePtr(complex double ***ptr, int imax, int jmax) {
  int i;
  for (i = 0; i < imax; i++) {
    freeComplexDoublePtr(ptr[i],jmax);
  }
  free(ptr);
  return;
}

void freeGrid(struct Grid *g) {
  // First, free all arrays:
  freeDoublePtr(ex, xSize);
  freeDoublePtr(ey, xSize + 1);
  freeDoublePtr(hz, xSize + 1);

  freeDoublePtr(heConst, xSize);
  freeDoublePtr(ehConst, xSize);
  freeDoublePtr(eqConst, xSize);
  freeDoublePtr(ABConst, xSize);
  freeComplexTriplePtr(qConst1, number_poles, xSize);
  freeComplexTriplePtr(qConst2, number_poles, xSize);
  freeComplexTriplePtr(qx, number_poles, xSize);
  freeComplexTriplePtr(qy, number_poles, xSize);
  freeDoublePtr(qxSum, xSize);
  freeDoublePtr(qySum, xSize);
  freeDoublePtr(iConst1, xSize);
  freeDoublePtr(iConst2, xSize);

  // Tridiagonal values:
  freeDoublePtr(a, xSize);
  freeDoublePtr(b, xSize);
  freeDoublePtr(c, xSize);

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
  free(reflDFT);
  free(tranDFT);
  free(qSumC);

  // And finally the grid itself:
  free(g);

  return;
}
