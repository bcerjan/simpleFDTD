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

/* Function for generating emscripten WASM file for embedding in a webpage. */

#include "fdtd_macro.h"
#include "fdtd_proto.h"
#include "array_proto.h"
#include "fdtd_field_proto.h"
#include "ezinc.h"
#include "sdl_funcs.h"
#include "emscripten.h"
#include <stdio.h>

// JS function call to update progress bar
extern void updateProgress(float percent);

// JS func to pass data to our chart:
extern void addReflData(double waveLength, double value);
extern void addTranData(double waveLength, double value);

// And to update the chart:
extern void updateChartData();

// To Enable "Run Simulation" button:
extern void enableButton();

void iterateSimulation(struct Grid *g) {
    // If first iteration, initialize SDL stuff:
    if (timeStep < 1) {
      imageInit(g);
    }
    //printf("Loop step: %i\timeStep",timeStep);
    HFieldUpdate(g, timeStep);
    EFieldUpdate(g);
    JFieldUpdate(g);
    lineSource(g, ABCSIZECONSTANT + 20, timeStep);
    //printf("ey at src: %f\n", ey[20][25]);
    DFTUpdate(g, timeStep);

    if (timeStep < 1500) {
      //printf("Loop step: %i\n",timeStep);
      //printf( "Drawing to screen...\n" );
      PlotField(g,2.5,0.0);
    }

    timeStep++;
    if (timeStep % 200 == 0) {
      updateProgress(100.0 * (float )timeStep / (float )maximumIteration);
    }

    // Check if we're done with the simulation:
    if (timeStep > maximumIteration) {
      // Finish the DFT calculations, accounting for empty run:
      finishFullDFT(g);

      for (int n = 0; n < NUMBERDFTFREQS; n++) {
        addReflData(wavelengthList[n]*dx*1e9, reflDFT[n]); // Convert back to nm from ppw
        addTranData(wavelengthList[n]*dx*1e9, tranDFT[n]); // ""
      }

      updateChartData();

      printf( "Finished loop\n" );
      updateProgress(100.0); // Finish Progress bar now that we're done
      freeGrid(g); // Free Grid struct
      imageFree(); // Free the SDL renderer
      enableButton();
      emscripten_cancel_main_loop(); // Cancel the loop when we run out of steps
    } /* ifCondition */
}

int fdtdSim(int metalChoice, int objectChoice, double objectSize, double environmentIndex, double objectIndex) {
  printf( "Started main...\n" );

  //struct Grid *g = malloc(sizeof(struct Grid));
  struct Grid *g;
  g = AllocateGridMemory();

  printf( "Allocated Grid\n" );

  InitializeFdtd(g, metalChoice, objectChoice, objectSize, environmentIndex, objectIndex); // First int for metal, second for object shape
  printf( "Initialized Grid\n" );
  // prepare matrix of object edges:
  findMatEdge(g);
  printf("Found Edges\n");

  // Initialize animation / SDL:


  maximumIteration = NUMBEROFITERATIONCONSTANT;
  interval = 0;
  timeStep = 0;

  emscripten_set_main_loop_arg(iterateSimulation, g, -1, 1); // call iterateSimulation, pass it data in g, let the browser determine framerate, and simulate an infinite loop

  return 0;
}
