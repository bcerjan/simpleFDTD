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

void freeTriplePtr(double ***ptr, int imax, int jmax) {
  int i;
  for (i = 0; i < imax; i++) {
    freeDoublePtr(ptr[i],jmax);
  }
  free(ptr);
  return;
}

void freeGrid(struct Grid *g) {
  // First, free all arrays:
  freeDoublePtr(ex, xSize);
  freeDoublePtr(ey, xSize + 1);
  freeDoublePtr(hz, xSize);
  freeDoublePtr(c1SumX, xSize);
  freeDoublePtr(c2SumX, xSize);
  freeDoublePtr(c1SumY, xSize);
  freeDoublePtr(c2SumY, xSize);
  freeDoublePtr(c3Sum, xSize);
  freeDoublePtr(c4Sum, xSize);
  freeDoublePtr(c5Sum, xSize);

  freeTriplePtr(c1Grid, number_poles, xSize);
  freeTriplePtr(c2Grid, number_poles, xSize);
  freeTriplePtr(c3Grid, number_poles, xSize);
  freeTriplePtr(c4Grid, number_poles, xSize);
  freeTriplePtr(c5Grid, number_poles, xSize);

  freeDoublePtr(dahz, xSize);
  freeDoublePtr(dbhz, xSize);

  freeTriplePtr(px, number_poles, xSize);
  freeTriplePtr(py, number_poles, xSize);
  freeTriplePtr(pxOld, number_poles, xSize);
  freeTriplePtr(pyOld, number_poles, xSize);

  freeDoublePtr(exOld, xSize);
  freeDoublePtr(eyOld, xSize + 1);
  freeDoublePtr(exOld2, xSize);
  freeDoublePtr(eyOld2, xSize + 1);
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
  free(hzOld);
  free(hGrad1);
  free(hGrad2);
  free(hGrad3);
  free(eGrad1);
  free(eGrad2);
  free(eGrad3);
  free(pmlSx);
  free(pmlSxOld);
  free(pmlSy);
  free(pmlSyOld);
  free(pmlTz);
  free(pmlTzOld);
  free(rx);
  free(rxOld);
  free(rxOld2);
  free(ry);
  free(ryOld);
  free(ryOld2);
  free(bz);
  free(bzOld);


  // And finally the grid itself:
  free(g);

  return;
}
