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

/* Function to run a simulation with a structure to test DFT results */

#include "fdtd_macro.h"
#include "fdtd_proto.h"
#include "array_proto.h"
#include "fdtd_field_proto.h"
#include "ezinc.h"
#include <stdio.h>

int main() {
  printf( "Started main...\n" );

  double exmax, exmin, eymax, eymin, jymax;
  double environmentIndex = 2.5;
  //struct Grid *g = malloc(sizeof(struct Grid));
  struct Grid *g;
  g = AllocateGridMemory();

  printf( "Allocated Grid\n" );

  InitializeFdtd(g, 0, -1, 20000.0, environmentIndex, 0.0); // First int for material, second for object shape, third for size, and fourth for dielectric environment
  printf( "Initialized Grid\n" );

  maximumIteration = NUMBEROFITERATIONCONSTANT;

  int n;
  int outInterval = 0;

  for (n = 0; n < maximumIteration; n++) {
//  for (n = 0; n < 15; n++) {
    HFieldUpdate(g, n);
    EFieldUpdate(g);
    JFieldUpdate(g);
    lineSource(g, ABCSIZECONSTANT + 20, n);
    //printf("ey at src: %f\n", ey[20][25]);
    DFTUpdate(g, n, environmentIndex);

    char tranEyFilename[100] = "test_output/structure_tran_raw_ey.h";
    FILE *tranEyDataPtr;

    // Write to header file for use later
    tranEyDataPtr = fopen(tranEyFilename, "a");
    fprintf(tranEyDataPtr, "%.17g,\n", ey[tranXPos][75]);
    fclose(tranEyDataPtr);

    char reflEyFilename[100] = "test_output/structure_refl_raw_ey.h";
    FILE *reflEyDataPtr;

    // Write to header file for use later
    reflEyDataPtr = fopen(reflEyFilename, "a");
    fprintf(reflEyDataPtr, "%.17g,\n", ey[reflXPos][75]);
    fclose(reflEyDataPtr);

    char tranHzFilename[100] = "test_output/structure_tran_raw_hz.h";
    FILE *tranHzDataPtr;

    // Write to header file for use later
    tranHzDataPtr = fopen(tranHzFilename, "a");
    fprintf(tranHzDataPtr, "%.17g,\n", hz[tranXPos][75]);
    fclose(tranHzDataPtr);

    char reflHzFilename[100] = "test_output/structure_refl_raw_hz.h";
    FILE *reflHzDataPtr;

    // Write to header file for use later
    reflHzDataPtr = fopen(reflHzFilename, "a");
    fprintf(reflHzDataPtr, "%.17g,\n", hz[reflXPos][75]);
    fclose(reflHzDataPtr);

    interval++;
  } /* nForLoop */

  // Scale our DFT's by the empty run:
  finishFullDFT(g);

  for (n = 0; n < NUMBERDFTFREQS; n++) {
    printf("reflDFT[%i]: %.17g\n",n,reflDFT[n] );
    //printf("tranDFT[%i]: %.17g\n",n,tranDFT[n] );
  }
  printf( "Finished loop\n" );

  freeGrid(g);
  return 0;
}
