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
#include <stdbool.h>
#include <math.h>

// JS function call to update progress bar
extern void updateProgress(float percent);

// JS func to pass data to our chart:
extern void addReflData(double waveLength, double value);
extern void addTranData(double waveLength, double value);

// And to update the chart:
extern void updateChartData();

// To Enable "Run Simulation" button:
extern void enableButton();

// Static variable to track if we're doing adaptive timing:
static bool adaptiveT = false;

void iterateSimulation(struct Grid *g) {
    // If first iteration, initialize SDL stuff and check for adaptive timing:
    if (timeStep < 1) {
      imageInit(g);
    } /* timeStep if Block */
    //printf("Loop step: %i\timeStep",timeStep);

    EFieldUpdate(g);
    lineSource(g, xSource, timeStep);
    QFieldUpdate(g);
    HFieldUpdate(g);


    DFTUpdate(g, timeStep);

    if (timeStep < 1500) {
      //printf("Loop step: %i\n",timeStep);
      //printf( "Drawing to screen...\n" );
      PlotField(g,2.5,0.0);
    }

    timeStep++;
    if (timeStep % 200 == 0) {
      if( adaptiveT > 0) {
        // This formula is a rough sigmoid to make the progress bar look a
        // little bit more active / more accurately represent how much time remains
        float percent = 1.0/(1.0 + (3000.0/(float ) timeStep)*expf(
                                  -((float )timeStep/3000.0 +
                                  1e-6/(float )AbsArrayMax(ey,xSize,ySize))
                                ) );
        updateProgress( 100.0*percent );
      } else {
        updateProgress(100.0 * (float )timeStep / (float )maximumIteration);
      } /* if/else block */
    }

    // Check for adaptive timing and if necessary increase maximumIteration
    // Only check every 50 time steps
    if (adaptiveT > 0 && timeStep % 50 == 0) {
      double maxEy = AbsArrayMax(ey,xSize,ySize);
      if (maxEy < 1e-6) { // Tolerance for field decay
        maximumIteration = 1; // This ends the loop on this iteration
      } else {
        maximumIteration += 51; // This indefinitely extends the loop
      } /* tolerance if block */
    } /* adaptiveT if block */

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

int fdtdSim(int metalChoice, int objectChoice, double objectXSize, double objectYSize,
  double environmentIndex, double objectIndex, bool adaptiveTime) {
  printf( "Started main...\n" );

  //struct Grid *g = malloc(sizeof(struct Grid));
  struct Grid *g;
  g = AllocateGridMemory();

  printf( "Allocated Grid\n" );

  InitializeFdtd(g, metalChoice, objectChoice, objectXSize, objectYSize,
    environmentIndex, objectIndex); // First int for metal, second for object shape

  printf( "Initialized Grid\n" );
  // prepare matrix of object edges:
  findMatEdge(g);
  printf("Found Edges\n");

  // Check if we're doing adaptive timing or not:
  if (adaptiveTime) {
    adaptiveT = true;
  } else {
    adaptiveT = false; // Need to reset it if we run a second time
  } /* ifCondition */

  maximumIteration = NUMBEROFITERATIONCONSTANT;

  interval = 0;
  timeStep = 0;

  emscripten_set_main_loop_arg(iterateSimulation, g, -1, 1); // call iterateSimulation, pass it data in g, let the browser determine framerate, and simulate an infinite loop

  return 0;
}
