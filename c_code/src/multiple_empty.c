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

/* Function to run empty simulations at a variety of background indices and
   write the output to a header file for use as normalization later on. */

#include "fdtd_macro.h"
#include "fdtd_proto.h"
#include "array_proto.h"
#include "fdtd_field_proto.h"
#include "ezinc.h"
#include "output_funcs.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
  printf( "Started main...\n" );

  double exmax, exmin, eymax, eymin, jymax;

  //struct Grid *g;

  int maximumIter = NUMBEROFITERATIONCONSTANT;

  int n,i,j;
  int outInterval = 0;
  double backInd, temp;
  // Divide by 10 to get actual background refractive index
  int minInd = 10;
  int maxInd = 12 + 1; // +1 is so you can write the max index you actually want
  //int minInd = 10;
  //int maxInd = 12;
  int numInd = maxInd - minInd;
  struct AuxIndexFields *Fields = malloc(sizeof(struct AuxIndexFields));

  // Arrays to store results as a function of background index
  Fields->nTranDFT = AllocateMemory(numInd, NUMBERDFTFREQS, 0.0);
  Fields->nReflDFT = AllocateMemory(numInd, NUMBERDFTFREQS, 0.0);
  Fields->nReEy = AllocateMemory(numInd, NUMBERDFTFREQS, 0.0);
  Fields->nImEy = AllocateMemory(numInd, NUMBERDFTFREQS, 0.0);
  Fields->nReHz = AllocateMemory(numInd, NUMBERDFTFREQS, 0.0);
  Fields->nImHz = AllocateMemory(numInd, NUMBERDFTFREQS, 0.0);

  // Number of indices
  Fields->numIndex = numInd;
  Fields->dftFreqs = NUMBERDFTFREQS;
  Fields->maxSteps = maximumIter;

  for (i = 0; i < numInd; i++ ) {
    struct Grid *g = malloc(sizeof(struct Grid));
    maximumIteration = NUMBEROFITERATIONCONSTANT;
    temp = (double  )(i + minInd);
    backInd = temp / 10.0;

    InitializeFdtd(g, 0, -1, 100.0, 100.0, backInd, 0.0); // First int for metal, second for object shape

    // Run simulation loop for this background index:
    for (n = 0; n < maximumIter; n++) {



      EFieldUpdate(g);
      lineSource(g, xSource, n);
      QFieldUpdate(g);
      HFieldUpdate(g);

      //printf("ey at src: %f\n", ey[20][25]);
      DFTUpdate(g, n);

      interval++;
    } /* nForLoop */

    // Finalize DFT
    finishEmptyDFT(g);
    printf("Writing to File...\n" );
    // Now store the values in our auxillary arrays:
    for (j = 0; j < NUMBERDFTFREQS; j++) {
      Fields->nTranDFT[i][j] = tranDFT[j];
      Fields->nReflDFT[i][j] = reflDFT[j];

      // As all of these are constant as a function of position, only store
      // them for one location:
      Fields->nReEy[i][j] = reEyReflDFT[j][ySize/2];
      Fields->nImEy[i][j] = imEyReflDFT[j][ySize/2];
      Fields->nReHz[i][j] = reHzReflDFT[j][ySize/2];
      Fields->nImHz[i][j] = imHzReflDFT[j][ySize/2];
    } /* jForLoop */

    printf( "Finished index loop: %i of %i\n", i, numInd - 1 );

    freeGrid(g);
  } /* iForLoop */

  WriteDFTFile(Fields);

  // Memory Deallocation:
  freeDoublePtr(Fields->nTranDFT, numInd);
  freeDoublePtr(Fields->nReflDFT, numInd);
  freeDoublePtr(Fields->nReEy, numInd);
  freeDoublePtr(Fields->nImEy, numInd);
  freeDoublePtr(Fields->nReHz, numInd);
  freeDoublePtr(Fields->nImHz, numInd);
  free(Fields);

  return 0;
}
