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

  //struct Grid *g = malloc(sizeof(struct Grid));
  struct Grid *g;
  g = AllocateGridMemory();

  printf( "Allocated Grid\n" );

  InitializeFdtd(g, 0, 0); // First int for metal, second for object shape
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
    DFTUpdate(g, n);

    char tranFilename[100] = "test_output/disk_tran_raw.h";
    FILE *tranDataPtr;

    // Write to header file for use later
    tranDataPtr = fopen(tranFilename, "a");
    fprintf(tranDataPtr, "%.17g,\n", ey[tranXPos][75]);
    fclose(tranDataPtr);

    char reflFilename[100] = "test_output/disk_refl_raw.h";
    FILE *reflDataPtr;

    // Write to header file for use later
    reflDataPtr = fopen(reflFilename, "a");
    fprintf(reflDataPtr, "%.17g,\n", ey[reflXPos][75]);
    fclose(reflDataPtr);

    interval++;
  } /* nForLoop */

  // Scale our DFT's by number of time steps:
  finishDFT(g);
  NormalizeDFT(g);

  for (n = 0; n < NUMBERDFTFREQS; n++) {
    printf("reflDFT[%i]: %.17g\n",n,reflDFT[n] );
    printf("tranDFT[%i]: %.17g\n",n,tranDFT[n] );
  }
  printf( "Finished loop\n" );
  //printf("%f\n", ey[70][10] );
  //WriteDFTFile(g);

  freeGrid(g);
  return 0;
}
