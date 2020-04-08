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

/* Function to run an empty simulation and write the output to a header file for use
   as normalization later on. */

#include "fdtd_macro.h"
#include "fdtd_proto.h"
#include "array_proto.h"
#include "fdtd_field_proto.h"
#include "ezinc.h"
#include <stdio.h>

int main() {
  printf( "Started main...\n" );

  double exmax, exmin, eymax, eymin, jymax;

  //struct Grid *g = malloc(sizeof(struct Grid));
  struct Grid *g;
  g = AllocateGridMemory();

  printf( "Allocated Grid\n" );

  InitializeFdtd(g, 0, -1, 100.0, 1.0, 0.0); // First int for metal, second for object shape
  printf( "Initialized Grid\n" );

  maximumIteration = NUMBEROFITERATIONCONSTANT;

  int n;
  int outInterval = 0;

  for (n = 0; n < maximumIteration; n++) {
  //for (n = 0; n < 15; n++) {
    HFieldUpdate(g, n);
    EFieldUpdate(g);
    JFieldUpdate(g);
    lineSource(g, ABCSIZECONSTANT + 20, n);
    //printf("ey at src: %f\n", ey[20][25]);
    DFTUpdate(g, n, 1.0);
    char tranFilename[100] = "test_output/empty_tran_raw.h";
    FILE *tranDataPtr;

    // Write to header file for use later
    tranDataPtr = fopen(tranFilename, "a");
    fprintf(tranDataPtr, "%.17g,\n", ey[tranXPos][75]);
    fclose(tranDataPtr);

    char reflFilename[100] = "test_output/empty_refl_raw.h";
    FILE *reflDataPtr;

    // Write to header file for use later
    reflDataPtr = fopen(reflFilename, "a");
    fprintf(reflDataPtr, "%.17g,\n", ey[reflXPos][75]);
    fclose(reflDataPtr);
    interval++;
  } /* nForLoop */

  // Finish taking the DFT:
  finishEmptyDFT(g);


  printf( "Finished loop\n" );
  // Output our fields:
  WriteDFTFile(g);

  freeGrid(g);
  return 0;
}
