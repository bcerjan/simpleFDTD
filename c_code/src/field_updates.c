#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "empty_refl_data.h"
#include "empty_tran_data.h"
#include "fdtd_macro.h"

void EFieldUpdate (struct Grid *g) {
  int i,j;
  double temporary;

  /***********************************************************************/
  //     Update electric fields (EX and EY)
  /***********************************************************************/
  // See Schneider 10.62 for equations
  for (i = 0; i < xSize; i++) {
    for (j = 1; j < ySize; j++) {        // j=0 = pec, so don't evaluate
      //ex[i][j] = caex[i][j] * ex[i][j] + cbex[i][j] * ( hz[i][j] - hz[i][j-1] ); /* Original update */
      exOld[i][j] = ex[i][j]; // Store previous field for polarization current
      temporary = (1.0/2.0) * (1.0 + cjj[i][j]) * jx[i][j];
      // If things go sideways, it might be the dx term on the H field... (same for y below)
      ex[i][j] = caex[i][j] * ex[i][j] + cbex[i][j] * ( (hz[i][j] - hz[i][j-1]) - temporary ); // Need dx on H fields?
    } /* jForLoop */
  } /* iForLoop */

  for (i = 1; i < xSize; i++) {            // i=0 = pec, so don't evaluate
    for (j = 0; j < ySize; j++) {
      //ey[i][j] = caey[i][j] * ey[i][j] + cbey[i][j] * ( hz[i-1][j] - hz[i][j] ); /* Original update */
      eyOld[i][j] = ey[i][j]; // Store previous field for polarization current
      temporary = (1.0/2.0) * (1.0 + cjj[i][j]) * jy[i][j];
      ey[i][j] = caey[i][j] * ey[i][j] + cbey[i][j] * ( (hz[i-1][j] - hz[i][j]) - temporary ); // Need dx on H fields?
    } /* jForLoop */
  } /* iForLoop */

  return;
}

void HFieldUpdate (struct Grid *g, int n) {
  int boundaryIndex,regionIndex,i,j,xStart,xStop,yStart,yStop;
  double hzx;

  /***********************************************************************/
  //     Update magnetic fields (HZ) in center (main) grid
  /***********************************************************************/


  regionIndex = 0;    // center (main) grid
  xStart = regionData[regionIndex].xStart;
  xStop  = regionData[regionIndex].xStop ;
  yStart = regionData[regionIndex].yStart;
  yStop  = regionData[regionIndex].yStop ;

  for (i = xStart; i < xStop; i++) {
    for (j = yStart; j < yStop; j++) {
      hz[i][j] = dahz[i][j] * hz[i][j] + dbhz[i][j] * ( ex[i][j+1] - ex[i][j] + ey[i][j] - ey[i+1][j] );
    } /* jForLoop */
  } /* iForLoop */


  /***********************************************************************/
  //     Update HZ in PML regions (hzx,hzy)
  /***********************************************************************/

  boundaryIndex = 0;
  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;
    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
        hzx = hz[i][j] - hzy[boundaryIndex];   // extract hzx
        hzx = dahz[i][j] * hzx + dbhz[i][j] * ( ey[i][j] - ey[i+1][j] );    // dahz,dbhz holds dahzx,dbhzx
        hzy[boundaryIndex] = dahzy[boundaryIndex] * hzy[boundaryIndex] + dbhzy[boundaryIndex] * ( ex[i][j+1] - ex[i][j] );
        hz[i][j] = hzx +  hzy[boundaryIndex];  // update hz
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  } /* forLoop */

  return;
}

void JFieldUpdate (struct Grid *g) { // I know, it's not actually a field, it's a current.
  int i,j;

 /***********************************************************************/
  //     Update polarization current (jx,jy) (actually are dx*jx or dx*jy)
  /***********************************************************************/
  for (i = 0; i < xSize; i++) {
    for (j = 1; j < ySize; j++) {        // j=0 = pec, so don't evaluate
      jx[i][j] = cjj[i][j] * dx * jx[i][j] + cje[i][j] * (ex[i][j] + exOld[i][j]); // need dx back?
    } /* jForLoop */
  } /* iForLoop */

  for (i = 1; i < xSize; i++) {            // i=0 = pec, so don't evaluate
    for (j = 0; j < ySize; j++) {
      jy[i][j] = cjj[i][j] * dx * jy[i][j] + cje[i][j] * (ey[i][j] + eyOld[i][j]); // need dx back?
    } /* jForLoop */
  } /* iForLoop */

  return;
}

void DFTUpdate (struct Grid *g, int n) {
  int regionIndex,yStart,yStop;
  int i,j;
  double temporary;
  double time, maxTime, kVal;
  // Need to convert ints to doubles
  maxTime = (double  )DFTPADDEDTIME;
  time = (double  )n;

  /***********************************************************************/
  //     Update DFT values
  /***********************************************************************/
  regionIndex = 0;    // center (main) grid
  yStart = regionData[regionIndex].yStart;
  yStop  = regionData[regionIndex].yStop ;

  /*char tranFilename[100] = "disk_tran_raw.h";
  FILE *tranDataPtr;

  // Write to header file for use later
  tranDataPtr = fopen(tranFilename, "a");
  fprintf(tranDataPtr, "%.17g,\n", ey[tranXPos][75]);
  fclose(tranDataPtr);*/

  for (i = 0; i < NUMBERDFTFREQS; i++) {
    for (j = yStart; j < yStop; j++) {
      kVal = (double )kList[i];
      //temporary = -1.0 * 0.5 * (ey[reflXPos][j]*hz[reflXPos][j]); // Poynting Flux through our line (positive x) (so we negate the result)
      temporary = ey[reflXPos][j];
      reReflDFT[i][j] += temporary * cos(2*pi*kVal*time/maxTime);
      imReflDFT[i][j] -= temporary * sin(2*pi*kVal*time/maxTime);

      temporary = ey[tranXPos][j];
      reTranDFT[i][j] += temporary * cos(2*pi*kVal*time/maxTime);
      imTranDFT[i][j] -= temporary * sin(2*pi*kVal*time/maxTime);
      //temporary = 0.5 * (ey[tranXPos][j]*hz[tranXPos][j]); // Poynting Flux through our line
    } /* jForLoop */
  } /* iForLoop */


  return;
}

// Function to write header files for our DFT normalization
void WriteDFTFile (struct Grid *g) {
  int i;
  char reflFilename[100] = "../include/fdtd/empty_refl_data.h";
  char tranFilename[100] = "../include/fdtd/empty_tran_data.h";
  FILE *reflDataPtr, *tranDataPtr;

  // Write to header file for use later
  reflDataPtr = fopen(reflFilename, "w");
  tranDataPtr = fopen(tranFilename, "w");

  // First add blocking definitions (just in case...)
  fprintf(reflDataPtr, "#ifndef REFL_EMPTY_DATA\n#define REFL_EMPTY_DATA\n");
  fprintf(tranDataPtr, "#ifndef TRAN_EMPTY_DATA\n#define TRAN_EMPTY_DATA\n");

  // Second add lines for number of freqs and time steps:
  fprintf(reflDataPtr, "#define reflSteps %i\n", maximumIteration);
  fprintf(reflDataPtr, "#define reflFreqs %i\n", NUMBERDFTFREQS);
  fprintf(tranDataPtr, "#define tranSteps %i\n", maximumIteration);
  fprintf(tranDataPtr, "#define tranFreqs %i\n", NUMBERDFTFREQS);

  // Now open our arrays:
  fprintf(reflDataPtr, "static const double emptyReflDFT[%i] = {\n", NUMBERDFTFREQS );
  fprintf(tranDataPtr, "static const double emptyTranDFT[%i] = {\n", NUMBERDFTFREQS );

  // Then define our arrays:

  for (i = 0; i < NUMBERDFTFREQS; i++) {
    fprintf(reflDataPtr, "%.17g,\n", reflDFT[i]);
    fprintf(tranDataPtr, "%.17g,\n", tranDFT[i]);
  } /* iForLoop */

  // End the arrays and if block:
  fprintf(reflDataPtr, "};\n#endif" );
  fprintf(tranDataPtr, "};\n#endif" );

  fclose(reflDataPtr);
  fclose(tranDataPtr);

  return;
}

// Function that normalizes DFT results based on stored data:
void NormalizeDFT (struct Grid *g) {
  int i;
  //double emptyReflDFT[NUMBERDFTFREQS];
  //double emptyTranDFT[NUMBERDFTFREQS];

  // Normalize our DFT's relative to empty run:
  for (i = 0; i < NUMBERDFTFREQS; i++) {
    reflDFT[i] = reflDFT[i] / ( emptyReflDFT[i] ) - 1.0; // -1.0 because there is always reflected light due to the source
    tranDFT[i] = tranDFT[i] / ( emptyTranDFT[i] );
  } /* iForLoop */

  return;
}

// Normalize by number of steps taken
void finishDFT (struct Grid *g) {
  int i,j;
  double temporary;

  int regionIndex = 0;    // center (main) grid
  int yStart = regionData[regionIndex].yStart;
  int yStop  = regionData[regionIndex].yStop ;

/*
  for (i = 0; i < NUMBERDFTFREQS; i++) {
    reflDFT[i] = reflDFT[i] / ( maximumIteration );
    tranDFT[i] = tranDFT[i] / ( maximumIteration );
  }*/ /* iForLoop */

  // For Direct flux calculation:
  // Apply scaling outside of the sum (Schneider 5.32-33)
  for (i = 0 ; i < NUMBERDFTFREQS; i++) {
    for (j = yStart; j < yStop; j++) {
      /*reEyDFT[i][j] = reEyDFT[i][j];
      imEyDFT[i][j] = imEyDFT[i][j];
      reHzDFT[i][j] = reHzDFT[i][j] / (maximumIteration);
      imHzDFT[i][j] = imHzDFT[i][j] / (maximumIteration - 1);*/

      // Now store the Poynting flux in positive x-direction as a function of position:
      reflDFT[i] += ( (reReflDFT[i][j] * reReflDFT[i][j]) + (imReflDFT[i][j] * imReflDFT[i][j]) );
      tranDFT[i] += ( (reTranDFT[i][j] * reTranDFT[i][j]) + (imTranDFT[i][j] * imTranDFT[i][j]) );

    } /* jForLoop */
  } /* iForLoop */

  return;
}
