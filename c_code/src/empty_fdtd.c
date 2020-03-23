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

  InitializeFdtd(g, 0, -1); // First int for metal, second for object shape
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
    /*if (outInterval == 50) {
      eymax = ArrayMax(ey, xSize, ySize);
      eymin = ArrayMin(ey, xSize, ySize);
      jymax = ArrayMax(jy, xSize, ySize);
      double time = n;
      printf("time step: %i\n", n);
      //printf("source value: %f\n", ezInc(time, 0.0));
      printf("eymax: %f\n", eymax);
      printf("eymin: %f\n", eymin);
      printf("jymax: %f\n", jymax);
      printf("jy[65][25]: %f\n", jy[65][25]);
      printf("-------------------\n");
      interval = 0;
    }*/
    interval++;
  } /* nForLoop */

  // Scale our DFT's by number of time steps:
  finishDFT(g);

  for (n = 0; n < NUMBERDFTFREQS; n++) {
    printf("reflDFT[%i]: %f\n",n,reflDFT[n] );
    printf("tranDFT[%i]: %f\n",n,tranDFT[n] );
  }
  printf( "Finished loop\n" );
  //printf("%f\n", ey[70][10] );
  WriteDFTFile(g);

  freeGrid(g);
  return 0;
}
