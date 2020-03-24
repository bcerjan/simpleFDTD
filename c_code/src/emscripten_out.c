/* Function for generating emscripten WASM file for embedding in a webpage. */

#include "fdtd_macro.h"
#include "fdtd_proto.h"
#include "array_proto.h"
#include "fdtd_field_proto.h"
#include "ezinc.h"
#include "sdl_funcs.h"
#include "emscripten.h"
#include <stdio.h>

void iterateSimulation(struct Grid *g) {
    //printf("Loop step: %i\timeStep",timeStep);
    HFieldUpdate(g, timeStep);
    EFieldUpdate(g);
    JFieldUpdate(g);
    lineSource(g, 30, timeStep);
    //printf("ey at src: %f\n", ey[20][25]);
    DFTUpdate(g, timeStep);

    if (timeStep < 350) {
      //printf("Loop step: %i\n",timeStep);
      //printf( "Drawing to screen...\n" );
      PlotField(g,3.0,0.0);
    }

    timeStep++;
    // Check if we're done with the simulation:
    if (timeStep > maximumIteration) {
      // Scale our DFT's by number of time steps:
      finishDFT(g);
      // And scale them based on an empty run:
      NormalizeDFT(g);

      for (int n = 0; n < NUMBERDFTFREQS; n++) {
        printf("reflDFT[%i]: %f\n",n,reflDFT[n] );
        printf("tranDFT[%i]: %f\n",n,tranDFT[n] );
      }
      printf( "Finished loop\n" );

      freeGrid(g);
      emscripten_cancel_main_loop(); // Cancel the loop when we run out of steps
    } /* ifCondition */
}


int fdtdSim(int metalChoice, int objectChoice) {
  printf( "Started main...\n" );

  //struct Grid *g = malloc(sizeof(struct Grid));
  struct Grid *g;
  g = AllocateGridMemory();

  printf( "Allocated Grid\n" );

  InitializeFdtd(g, metalChoice, objectChoice); // First int for metal, second for object shape
  printf( "Initialized Grid\n" );
  // prepare matrix of object edges:
  findMatEdge(g);
  printf("Found Edges\n");
  //printf("max edgeMat: %f\n", ArrayMax(edgeMat,xSize,ySize));
  //printf("min edgeMat: %f\n", ArrayMin(edgeMat,xSize,ySize));
  maximumIteration = NUMBEROFITERATIONCONSTANT;
  interval = 0;
  timeStep = 0;

  emscripten_set_main_loop_arg(iterateSimulation, g, -1, 1); // call iterateSimulation, pass it data in g, let the browser determine framerate, and simulate an infinite loop

  return 0;
}
