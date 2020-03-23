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
  freeDoublePtr(object_locs, xSize);

  freeDoublePtr(reReflDFT, NUMBERDFTFREQS);
  freeDoublePtr(imReflDFT, NUMBERDFTFREQS);
  freeDoublePtr(reTranDFT, NUMBERDFTFREQS);
  freeDoublePtr(imTranDFT, NUMBERDFTFREQS);

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
