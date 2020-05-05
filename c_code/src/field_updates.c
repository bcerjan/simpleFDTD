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
#include "fdtd_proto.h"

void EFieldUpdate (struct Grid *g) {
  int i,j,p;

  // See Prokopidis and Zografopoulos eq. 23 with Pd and all D terms set to 0.
  for (i = 0; i < xSize; i++) {
    for (j = 1; j < ySize; j++) {        // j=0 = pec, so don't evaluate
      exOld2[i][j] = exOld[i][j]; // E at n - 1
      exOld[i][j] = ex[i][j]; // Store previous field for polarization

      ex[i][j] = ( dt * (hz[i][j] - hz[i][j-1]) + \
        c4Sum[i][j] * ex[i][j] - c5Sum[i][j] * exOld2[i][j] - \
        c1SumX[i][j] - c2SumX[i][j] ) / c3Sum[i][j];
    } /* jForLoop */
  } /* iForLoop */

  for (i = 1; i < xSize; i++) {            // i=0 = pec, so don't evaluate
    for (j = 0; j < ySize; j++) {
      eyOld2[i][j] = eyOld[i][j];
      eyOld[i][j] = ey[i][j]; // Store previous field for polarization current

      ey[i][j] = ( dt * (hz[i][j] - hz[i][j-1]) + \
        c4Sum[i][j] * ey[i][j] - c5Sum[i][j] * eyOld2[i][j] - \
        c1SumY[i][j] - c2SumY[i][j] ) / c3Sum[i][j];
    } /* jForLoop */
  } /* iForLoop */

  return;
}

void HFieldUpdate (struct Grid *g, int n) {
  int i,j;

  /***********************************************************************/
  //     Update magnetic fields (HZ) in center (main) grid
  /***********************************************************************/

  for (i = 0; i < xSize - 1; i++) {
    for (j = 0; j < ySize - 1; j++) {
      hz[i][j] = dahz[i][j] * hz[i][j] + dbhz[i][j] * ( ex[i][j+1] - ex[i][j] + ey[i][j] - ey[i+1][j] );
    } /* jForLoop */
  } /* iForLoop */

  return;
}

// Auxiliary fields for PML (both E and H)
// Call after E-field/H-Field update
// Currently using E and H as stand-ins for R and B in the PML region as they
// are getting updated in this region the same way R and B would be.
void SFieldUpdate (struct Grid *g) {
  int i,j,p,regionIndex,boundaryIndex,xStop,xStart,yStop,yStart;

  boundaryIndex = 0;
  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;
    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
        // Store previous versions for E / H updates
        pmlSxOld[boundaryIndex] = pmlSx[boundaryIndex];
        pmlSyOld[boundaryIndex] = pmlSy[boundaryIndex];
        pmlTzOld[boundaryIndex] = pmlTz[boundaryIndex];
        pmlSx[boundaryIndex] = eGrad1[boundaryIndex]*pmlSx[boundaryIndex] + \
                               eGrad2[boundaryIndex]*rx[boundaryIndex] - \
                               eGrad3[boundaryIndex]*rxOld[boundaryIndex];
        pmlSy[boundaryIndex] = eGrad1[boundaryIndex]*pmlSy[boundaryIndex] + \
                               eGrad2[boundaryIndex]*ry[boundaryIndex] - \
                               eGrad3[boundaryIndex]*ryOld[boundaryIndex];
        pmlTz[boundaryIndex] = hGrad1[boundaryIndex]*pmlTz[boundaryIndex] + \
                               hGrad2[boundaryIndex]*bz[boundaryIndex] - \
                               hGrad3[boundaryIndex]*bzOld[boundaryIndex];

        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  } /* region forLoop */
  return;
}

// Update E and H in PML regions (equivalently, find R and B)
// Call after SFieldUpdate
void PMLFieldUpdate (struct Grid *g) {
  int i,j,p,regionIndex,boundaryIndex,xStop,xStart,yStop,yStart;

  boundaryIndex = 0;
  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;
    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
        ex[i][j] = eGrad1[boundaryIndex]*exOld[i][j] + \
                   eGrad2[boundaryIndex]*pmlSx[boundaryIndex] - \
                   eGrad3[boundaryIndex]*pmlSxOld[boundaryIndex];
        ey[i][j] = eGrad1[boundaryIndex]*eyOld[i][j] + \
                   eGrad2[boundaryIndex]*pmlSy[boundaryIndex] - \
                   eGrad3[boundaryIndex]*pmlSyOld[boundaryIndex];
        hz[i][j] = hGrad1[boundaryIndex]*hzOld[boundaryIndex] + \
                   hGrad2[boundaryIndex]*pmlTz[boundaryIndex] - \
                   hGrad3[boundaryIndex]*pmlTzOld[boundaryIndex];
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  } /* region forLoop */
  return;
}

// One of several auxiliary fields for the PML:
void RFieldUpdate (struct Grid *g) {
  int i,j,p,regionIndex,boundaryIndex;

  boundaryIndex = 0;
  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;
    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
      rxOld2[boundaryIndex] = rxOld[boundaryIndex]; // E at n - 2
      rxOld[boundaryIndex] = rx[boundaryIndex]; // Store previous field for polarization

      rx[boundaryIndex] = ( dt * (hz[i][j+1] - hz[i][j]) + \
        c4Sum[i][j] * rx[boundaryIndex] - c5Sum[i][j] * rxOld2[boundaryIndex] - \
        c1SumX[i][j] - c2SumX[i][j] ) / c3Sum[i][j];

      ryOld2[boundaryIndex] = ryOld[boundaryIndex];
      ryOld[boundaryIndex] = ry[boundaryIndex]; // Store previous field for polarization current

      ry[boundaryIndex] = ( dt * (hz[i][j+1] - hz[i][j]) + \
        c4Sum[i][j] * ry[boundaryIndex] - c5Sum[i][j] * ryOld2[boundaryIndex] - \
        c1SumY[i][j] - c2SumY[i][j] ) / c3Sum[i][j];
      boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  }
  return;
}

void BFieldUpdate (struct Grid *g, int n) {
  int i,j,regionIndex,boundaryIndex;

  /***********************************************************************/
  //     Update magnetic fields (HZ) in center (main) grid
  /***********************************************************************/

  boundaryIndex = 0;
  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;
    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
      bzOld[boundaryIndex] = bz[boundaryIndex];
      bz[boundaryIndex] = dahz[i][j] * bz[boundaryIndex] + dbhz[i][j] * ( ex[i][j+1] - ex[i][j] + ey[i][j] - ey[i+1][j] );
      boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  } /* region forLoop */

  return;
}

// This should be called _after_ EFieldUpdate
void PFieldUpdate (struct Grid *g) { // I know, it's not actually a field, it's the polarization.
  int i,j,p;
  double tempOld;

  for (p = 0; p < number_poles; p++) {
    for (i = 0; i < xSize; i++) {
      for (j = 1; j < ySize; j++) { // j=0 -> pec
        tempOld = pxOld[p][i][j];
        pxOld[p][i][j] = px[p][i][j];
        px[p][i][j] = c1Grid[p][i][j]*px[p][i][j] + \
         c2Grid[p][i][j]*tempOld + c3Grid[p][i][j]*ex[i][j] + \
         c4Grid[p][i][j]*exOld[i][j] + c5Grid[p][i][j]*exOld2[i][j];

        if( p < 1 ) { // since we need to reset from previous time steps
          c1SumX[i][j] = (1.0 - c1Grid[p][i][j])*px[p][i][j];
          c2SumX[i][j] = c2Grid[p][i][j]*pxOld[p][i][j];
        } else { // Now we're doing our running sum:
          c1SumX[i][j] += (1.0 - c1Grid[p][i][j])*px[p][i][j];
          c2SumX[i][j] += c2Grid[p][i][j]*pxOld[p][i][j];
        } /* ifBlock */
      } /* jForLoop */
    } /* iForLoop */
  } /* pForLoop */

  for (p = 0; p < number_poles; p++) {
    for (i = 1; i < xSize; i++) {            // i=0 = pec, so don't evaluate
      for (j = 0; j < ySize; j++) {
        tempOld = pyOld[p][i][j];
        pyOld[p][i][j] = py[p][i][j];
        py[p][i][j] = c1Grid[p][i][j]*py[p][i][j] + \
         c2Grid[p][i][j]*tempOld + c3Grid[p][i][j]*ey[i][j] + \
         c4Grid[p][i][j]*eyOld[i][j] + c5Grid[p][i][j]*eyOld2[i][j];

         if( p < 1 ) { // since we need to reset from previous time steps
           c1SumY[i][j] = (1.0 - c1Grid[p][i][j])*py[p][i][j];
           c2SumY[i][j] = c2Grid[p][i][j]*pyOld[p][i][j];
         } else { // Now we're doing our running sum:
           c1SumY[i][j] += (1.0 - c1Grid[p][i][j])*py[p][i][j];
           c2SumY[i][j] += c2Grid[p][i][j]*pyOld[p][i][j];
         } /* ifBlock */
      } /* jForLoop */
    } /* iForLoop */
  } /* pForLoop */
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
