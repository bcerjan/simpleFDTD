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
#include "output_funcs.h"
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


  printf("AbsMax c3Sum: %.17g\n", AbsArrayMax(c3Sum,xSize,ySize));
  printf("AbsMax c4Sum: %.17g\n", AbsArrayMax(c4Sum,xSize,ySize));
  printf("AbsMax c5Sum: %.17g\n", AbsArrayMax(c5Sum,xSize,ySize));
  printf("---------------------\n");

  maximumIteration = NUMBEROFITERATIONCONSTANT;

  int n;
  int outInterval = 0;

  //for (n = 0; n < maximumIteration; n++) {
  for (n = 0; n < 31; n++) {

    HFieldUpdate(g);
    //printf( "H Updated\n" );
    EFieldUpdate(g);
    //printf( "E Updated\n" );
    BFieldUpdate(g);
    //printf( "B Updated\n" );
    PFieldUpdate(g);
    //printf( "P Updated\n" );
    RFieldUpdate(g);
    //printf( "R Updated\n" );
    SFieldUpdate(g);
    //printf( "S Updated\n" );
    PMLFieldUpdate(g);
    //printf( "PML Updated\n" );
    lineSource(g, xSource, n);
    //printf("ey at src: %f\n", ey[20][25]);
    DFTUpdate(g, n);
    if( n % 5 == 0 ){
    printf("n: %i\n",n);
    printf("AbsMax py: %.17g\n", AbsArrayMax(py[0],xSize,ySize));
    printf("AbsMax ry: %.17g\n", AbsVectorMax(ry,9056));
    printf("AbsMax bz: %.17g\n", AbsVectorMax(bz,9056));
    printf("AbsMax c1SumY: %.17g\n", AbsArrayMax(c1SumY,xSize,ySize));
    printf("AbsMax c2SumY: %.17g\n", AbsArrayMax(c2SumY,xSize,ySize));
    printf("AbsMax ey: %.17g\n", AbsArrayMax(ey,xSize,ySize));
    //printf("AbsMax sy: %.17g\n", AbsVectorMax(pmlSy,9056));
    //printf("AbsMax tz: %.17g\n", AbsVectorMax(pmlTz,9056));
    printf("---------------------\n");}
    /*char tranEyFilename[100] = "test_output/empty_tran_raw_ey.h";
    FILE *tranEyDataPtr;

    // Write to header file for use later
    tranEyDataPtr = fopen(tranEyFilename, "a");
    fprintf(tranEyDataPtr, "%.17g,\n", ey[tranXPos][75]);
    fclose(tranEyDataPtr);

    char reflEyFilename[100] = "test_output/empty_refl_raw_ey.h";
    FILE *reflEyDataPtr;

    // Write to header file for use later
    reflEyDataPtr = fopen(reflEyFilename, "a");
    fprintf(reflEyDataPtr, "%.17g,\n", ey[reflXPos][75]);
    fclose(reflEyDataPtr);

    char tranHzFilename[100] = "test_output/empty_tran_raw_hz.h";
    FILE *tranHzDataPtr;

    // Write to header file for use later
    tranHzDataPtr = fopen(tranHzFilename, "a");
    fprintf(tranHzDataPtr, "%.17g,\n", hz[tranXPos][75]);
    fclose(tranHzDataPtr);

    char reflHzFilename[100] = "test_output/empty_refl_raw_hz.h";
    FILE *reflHzDataPtr;

    // Write to header file for use later
    reflHzDataPtr = fopen(reflHzFilename, "a");
    fprintf(reflHzDataPtr, "%.17g,\n", hz[reflXPos][75]);
    fclose(reflHzDataPtr);*/

    interval++;
  } /* nForLoop */

  // Finish taking the DFT:
  finishEmptyDFT(g);


  printf( "Finished loop\n" );
  // Output our fields:
  //WriteDFTFile(g);

  freeGrid(g);
  return 0;
}
