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

  for (i = 0; i < NUMBERDFTFREQS; i++) {
    for (j = yStart; j < yStop; j++) {
      kVal = (double )kList[i];

      // First do Ey accumulations:
      temporary = ey[reflXPos][j];
      reEyReflDFT[i][j] += temporary * cos(2*pi*kVal*time/maxTime);
      imEyReflDFT[i][j] -= temporary * sin(2*pi*kVal*time/maxTime);

      temporary = ey[tranXPos][j];
      reEyTranDFT[i][j] += temporary * cos(2*pi*kVal*time/maxTime);
      imEyTranDFT[i][j] -= temporary * sin(2*pi*kVal*time/maxTime);

      // Then Hz accumulations:
      temporary = hz[reflXPos][j];
      reHzReflDFT[i][j] += temporary * cos(2*pi*kVal*time/maxTime);
      imHzReflDFT[i][j] -= temporary * sin(2*pi*kVal*time/maxTime);

      temporary = hz[tranXPos][j];
      reHzTranDFT[i][j] += temporary * cos(2*pi*kVal*time/maxTime);
      imHzTranDFT[i][j] -= temporary * sin(2*pi*kVal*time/maxTime);

    } /* jForLoop */
  } /* iForLoop */


  return;
}

// Function to flatten the phase profile of an input FFT spectrum
// Only currently used in finishEmptyDFT and finishFullDFT
// We could honestly replace this with just complexField[i][j] == abs(complexField[i][j])
// because that's what we're doing anyway by setting the phase to 0.
// ----> Is this a problem? <-----
// Not used right now, as it breaks things (and is probably a bad idea anyway...)
void flattenPhase(double **reField, double **imField, int numFreqs,
  int length, double targetPhase) {
  int i,j;
  double phase,reComp,imComp,phaseAdjustment;

  for (i = 0; i < numFreqs; i++) {
    for (j = 0; j < length; j++) {
      // Fix current components so we can update actual fields in place
      reComp = reField[i][j];
      imComp = imField[i][j];
      // Calculate phase angle at this frequency and location
      phase = atan2(imComp, reComp);
      phaseAdjustment = targetPhase - phase;
      // Flatten phase profile so we can use a single calibration run
      // These are expansions of: F[f,pos] * exp( -i * phase ) <- note the minus
      reField[i][j] = reComp * cos(phaseAdjustment) + imComp * sin(phaseAdjustment);
      imField[i][j] = imComp * cos(phaseAdjustment) - reComp * sin(phaseAdjustment);
    } /* jForLoop */
  } /* iForLoop */
  return;
}

// Function that normalizes DFT results based on stored data for both
// reflected and transmitted fields:
void finishFullDFT (struct Grid *g) {
  int i,j;
  int regionIndex = 0;    // center (main) grid
  int yStart = regionData[regionIndex].yStart;
  int yStop  = regionData[regionIndex].yStop ;
  double reEy,imEy,reHz,imHz;
  /**double emptyReEyRefl[NUMBERDFTFREQS][ySize],\
  emptyImEyRefl[NUMBERDFTFREQS][ySize],\
  emptyReHzRefl[NUMBERDFTFREQS][ySize],\
  emptyImHyRefl[NUMBERDFTFREQS][ySize];**/

  // First do Transmission as it is simpler:
  for (i = 0; i < NUMBERDFTFREQS; i++) {
    for (j = yStart; j < yStop; j++) {
      tranDFT[i] += ( (reEyTranDFT[i][j] * reHzTranDFT[i][j]) + (imHzTranDFT[i][j] * imEyTranDFT[i][j]) );
    } /* jForLoop */

    /* Now we normalize based on the empty run data for the background requested. */
    tranDFT[i] = tranDFT[i] / ( emptyTranDFT[refractiveIndexIndex][i] );
  } /* iForLoop */

  // Now, we need to compute the reflected flux accounting for the intial fields:
  // See: https://meep.readthedocs.io/en/latest/Introduction/#transmittancereflectance-spectra

  for (i = 0; i < NUMBERDFTFREQS; i++) {
    for (j = yStart; j < yStop; j++) {
      reEy = reEyReflDFT[i][j] - emptyReEyRefl[refractiveIndexIndex][i];
      imEy = imEyReflDFT[i][j] - emptyImEyRefl[refractiveIndexIndex][i];
      reHz = reHzReflDFT[i][j] - emptyReHzRefl[refractiveIndexIndex][i];
      imHz = imHzReflDFT[i][j] - emptyImHzRefl[refractiveIndexIndex][i];
      reflDFT[i] -=  ( (reEy * reHz) + (imEy * imHz) ); // -= is because we are pointing in the negative x direction
    } /* jForLoop */
    reflDFT[i] = reflDFT[i] / emptyReflDFT[refractiveIndexIndex][i]; // Normalize to empty run
  } /* iForLoop */

  return;
}

// Convert Ey(w) and Hz(w) to P(w) for transmitted fields:
void finishEmptyDFT (struct Grid *g) {
  int i,j;

  int regionIndex = 0;    // center (main) grid
  int yStart = regionData[regionIndex].yStart;
  int yStop  = regionData[regionIndex].yStop ;

  // Poynting flux calculation Integral of (Ey* x Hz):
  for (i = 0 ; i < NUMBERDFTFREQS; i++) {
    for (j = yStart; j < yStop; j++) {
      tranDFT[i] += ( (reEyTranDFT[i][j] * reHzTranDFT[i][j]) + (imHzTranDFT[i][j] * imEyTranDFT[i][j]) );
      reflDFT[i] -= ( (reEyReflDFT[i][j] * reHzReflDFT[i][j]) + (imHzReflDFT[i][j] * imEyReflDFT[i][j]) ); // Minus as we want -X direction
    } /* jForLoop */
  } /* iForLoop */

  return;
}
