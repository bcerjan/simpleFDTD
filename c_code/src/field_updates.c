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


/* Function to store R, S, and P fields at the top of each timestep: */
void StoreFields (struct Grid *g) {
  int i,j,p,regionIndex,boundaryIndex,xStop,xStart,yStop,yStart;

  boundaryIndex = 0;
  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;
    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
        /* S/T updates (maybe T is worthless?) */
        pmlSxOld[boundaryIndex] = pmlSx[boundaryIndex]; // These are all now at time step "n"
        pmlSyOld[boundaryIndex] = pmlSy[boundaryIndex];
        pmlTzOld[boundaryIndex] = pmlTz[boundaryIndex];

        /* R storage */
        rxOld2[boundaryIndex] = rxOld[boundaryIndex]; // Rx at n - 1
        rxOld[boundaryIndex] = rx[boundaryIndex]; // Rx at n
        ryOld2[boundaryIndex] = ryOld[boundaryIndex]; // Ry at n - 1
        ryOld[boundaryIndex] = ry[boundaryIndex]; // Ry at n
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  } /* region forLoop */

  /* Store Px / Py */
  for (p = 0; p < number_poles; p++) {
    for (i = 0; i < xSize; i++) {
      for (j = 0; j < ySize; j++) { // j=0 -> pec
        pxOld2[p][i][j] = pxOld[p][i][j]; // Px at n - 1
        pxOld[p][i][j] = px[p][i][j]; // Px at n
        pyOld2[p][i][j] = pyOld[p][i][j]; // Py at n - 1
        pyOld[p][i][j] = py[p][i][j]; // Py at n
      } /* jForLoop */
    } /* iForLoop */
  } /* pForLoop */

  return;
}


// Auxiliary "E" field in PML region:
void RFieldUpdate (struct Grid *g) {
  int i,j,regionIndex,boundaryIndex,xStop,xStart,yStop,yStart;
  double temphz1,temphz2;

  boundaryIndex = 0;
  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;

    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
        if ( j < 1 ) { // To avoid trying to access out of bounds memory
          temphz1 = hz[i][j];
          //temphz1 = 0.0;
        } else {
          temphz1 = hz[i][j-1];
        } /* ifBlock */

        if ( i < 1 ) {
          temphz2 = hz[i][j];
          //temphz2 = 0.0;
        } else {
          temphz2 = hz[i-1][j];
        }
        /*if ( i < 1 && j < 1 ) {
          printf("H Term: %.17g\n", dt*(temphz2 - hz[i][j])/c3Sum[i][j]);
          printf("c4 Term: %.17g\n", c4Sum[i][j]*rx[boundaryIndex]/c3Sum[i][j]);
          printf("c5 Term: %.17g\n", -1.0*c5Sum[i][j]*rxOld2[boundaryIndex]/c3Sum[i][j]);
          printf("c1 Term: %.17g\n", c1SumY[i][j]/c3Sum[i][j]);
          printf("c2 Term: %.17g\n", c2SumY[i][j]/c3Sum[i][j]);
          printf("Numerator: %.17g\n", c4Sum[i][j]*rx[boundaryIndex] - c5Sum[i][j]*rxOld2[boundaryIndex] - c1SumY[i][j] - c2SumY[i][j] );
        }*/
        rx[boundaryIndex] = ( (dt/dx) * (hz[i][j] - temphz1) + \
          c4Sum[i][j] * rx[boundaryIndex] - c5Sum[i][j] * rxOld2[boundaryIndex] - \
          c1SumX[i][j] - c2SumX[i][j] ) / c3Sum[i][j];

        ry[boundaryIndex] = ( (dt/dx) * (temphz2 - hz[i][j]) + \
          c4Sum[i][j] * ry[boundaryIndex] - c5Sum[i][j] * ryOld2[boundaryIndex] - \
          c1SumY[i][j] - c2SumY[i][j] ) / c3Sum[i][j];
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  } /* region forLoop */
  return;
}

// Update Px/Py in PML region:
void PPMLFieldUpdate (struct Grid *g) { // I know, it's not actually a field, it's the polarization.
  int i,j,p,regionIndex,boundaryIndex,xStart,xStop,yStart,yStop;

  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;
    /* Px */
    for (p = 0; p < number_poles; p++) {
      boundaryIndex = 0; // need to reset our counter for each pole
      for (i = xStart; i < xStop; i++) {
        for (j = yStart; j < yStop; j++) { // i=0 = pec, so don't evaluate
          px[p][i][j] = c1Grid[p][i][j]*px[p][i][j] + \
           c2Grid[p][i][j]*pxOld2[p][i][j] + c3Grid[p][i][j]*rx[boundaryIndex] + \
           c4Grid[p][i][j]*rxOld[boundaryIndex] + c5Grid[p][i][j]*rxOld2[boundaryIndex];

          py[p][i][j] = c1Grid[p][i][j]*py[p][i][j] + \
           c2Grid[p][i][j]*pyOld2[p][i][j] + c3Grid[p][i][j]*ry[boundaryIndex] + \
           c4Grid[p][i][j]*ryOld[boundaryIndex] + c5Grid[p][i][j]*ryOld2[boundaryIndex];

          if( p < 1 ) { // since we need to reset from previous time steps
            c1SumX[i][j] = (c1Grid[p][i][j] - 1.0)*px[p][i][j];
            c2SumX[i][j] = c2Grid[p][i][j]*pxOld2[p][i][j];
            c1SumY[i][j] = (c1Grid[p][i][j] - 1.0)*py[p][i][j];
            c2SumY[i][j] = c2Grid[p][i][j]*pyOld2[p][i][j];
          } else { // Now we're doing our running sum:
            c1SumX[i][j] += (c1Grid[p][i][j] - 1.0)*px[p][i][j];
            c2SumX[i][j] += c2Grid[p][i][j]*pxOld2[p][i][j];
            c1SumY[i][j] += (c1Grid[p][i][j] - 1.0)*py[p][i][j];
            c2SumY[i][j] += c2Grid[p][i][j]*pyOld2[p][i][j];
          } /* ifBlock */

          boundaryIndex++;
        } /* jForLoop */
      } /* iForLoop */
    } /* pForLoop */
  } /* region forLoop */
  return;
}

// Auxiliary field for updating E in PML region
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
        pmlSx[boundaryIndex] = sGrad1[boundaryIndex]*pmlSx[boundaryIndex] + \
                               sGrad2[boundaryIndex]*rx[boundaryIndex] - \
                               sGrad3[boundaryIndex]*rxOld[boundaryIndex];
        pmlSy[boundaryIndex] = sGrad1[boundaryIndex]*pmlSy[boundaryIndex] + \
                               sGrad2[boundaryIndex]*ry[boundaryIndex] - \
                               sGrad3[boundaryIndex]*ryOld[boundaryIndex];
        /*pmlTz[boundaryIndex] = tGrad1[boundaryIndex]*pmlTz[boundaryIndex] + \
                               tGrad2[boundaryIndex]*bz[boundaryIndex] - \
                               tGrad3[boundaryIndex]*bzOld[boundaryIndex];*/
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  } /* region forLoop */
  return;
}


void EFieldUpdate (struct Grid *g) {
  int i,j,xStart,xStop,yStart,yStop,boundaryIndex,regionIndex;


  /* Update E in main region */
  xStart = regionData[0].xStart;
  xStop  = regionData[0].xStop ;
  yStart = regionData[0].yStart;
  yStop  = regionData[0].yStop ;

  // See Prokopidis and Zografopoulos eq. 23 with Pd and all D terms set to 0.
  for (i = xStart; i < xStop; i++) {
    for (j = yStart; j < yStop; j++) {        // j=0 = pec, so don't evaluate
      exOld2[i][j] = exOld[i][j]; // E at n - 1
      exOld[i][j] = ex[i][j]; // Store previous field for polarization (E at n)

      ex[i][j] = ( (dt/dx) * (hz[i][j] - hz[i][j-1]) + \
        c4Sum[i][j] * ex[i][j] - c5Sum[i][j] * exOld2[i][j] - \
        c1SumX[i][j] - c2SumX[i][j] ) / c3Sum[i][j]; // E at n + 1

      eyOld2[i][j] = eyOld[i][j];
      eyOld[i][j] = ey[i][j]; // Store previous field for polarization current

      ey[i][j] = ( (dt/dx) * (hz[i-1][j] - hz[i][j]) + \
        c4Sum[i][j] * ey[i][j] - c5Sum[i][j] * eyOld2[i][j] - \
        c1SumY[i][j] - c2SumY[i][j] ) / c3Sum[i][j];
    } /* jForLoop */
  } /* iForLoop */

  /* Update E in PML */
  boundaryIndex = 0;
  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;
    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
        ex[i][j] = eGrad1[boundaryIndex]*ex[i][j] + \
                   eGrad2[boundaryIndex]*pmlSx[boundaryIndex] - \
                   eGrad3[boundaryIndex]*pmlSxOld[boundaryIndex];
        ey[i][j] = eGrad1[boundaryIndex]*ey[i][j] + \
                   eGrad2[boundaryIndex]*pmlSy[boundaryIndex] - \
                   eGrad3[boundaryIndex]*pmlSyOld[boundaryIndex];
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  } /* region forLoop */

  return;
}

// Update Px/Py in main grid:
void PFieldUpdate (struct Grid *g) { // I know, it's not actually a field, it's the polarization.
  int i,j,p,xStop,xStart,yStop,yStart;
  xStart = regionData[0].xStart;
  xStop  = regionData[0].xStop ;
  yStart = regionData[0].yStart;
  yStop  = regionData[0].yStop ;

  for (p = 0; p < number_poles; p++) {
    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) { // i=0 = pec, so don't evaluate
        px[p][i][j] = c1Grid[p][i][j]*px[p][i][j] + \
         c2Grid[p][i][j]*pxOld2[p][i][j] + c3Grid[p][i][j]*ex[i][j] + \
         c4Grid[p][i][j]*exOld[i][j] + c5Grid[p][i][j]*exOld2[i][j];

         py[p][i][j] = c1Grid[p][i][j]*py[p][i][j] + \
          c2Grid[p][i][j]*pyOld2[p][i][j] + c3Grid[p][i][j]*ey[i][j] + \
          c4Grid[p][i][j]*eyOld[i][j] + c5Grid[p][i][j]*eyOld2[i][j];

         if( p < 1 ) { // since we need to reset from previous time steps
           c1SumX[i][j] = (c1Grid[p][i][j] - 1.0)*px[p][i][j];
           c2SumX[i][j] = c2Grid[p][i][j]*pxOld2[p][i][j];
           c1SumY[i][j] = (c1Grid[p][i][j] - 1.0)*py[p][i][j];
           c2SumY[i][j] = c2Grid[p][i][j]*pyOld2[p][i][j];
         } else { // Now we're doing our running sum:
           c1SumX[i][j] += (c1Grid[p][i][j] - 1.0)*px[p][i][j];
           c2SumX[i][j] += c2Grid[p][i][j]*pxOld2[p][i][j];
           c1SumY[i][j] += (c1Grid[p][i][j] - 1.0)*py[p][i][j];
           c2SumY[i][j] += c2Grid[p][i][j]*pyOld2[p][i][j];
         } /* ifBlock */
      } /* jForLoop */
    } /* iForLoop */
  } /* pForLoop */
  return;
}

// Update stand-in for H in PML (B)
void BFieldUpdate (struct Grid *g) {
  int i,j,regionIndex,boundaryIndex,xStop,xStart,yStop,yStart;

  boundaryIndex = 0;
  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;
    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
      bzOld2[boundaryIndex] = bzOld[boundaryIndex]; // Bz at n - 1
      bzOld[boundaryIndex] = bz[boundaryIndex]; // Bz at n
      bz[boundaryIndex] = -1.0*bGrad1[boundaryIndex] * ( ex[i][j+1] - ex[i][j] + ey[i][j] - ey[i+1][j] ) - \
        bGrad2[boundaryIndex] * ( exOld[i][j+1] - exOld[i][j] + eyOld[i][j] - eyOld[i+1][j] ) + \
        bGrad3[boundaryIndex] * bz[boundaryIndex] - bGrad4[boundaryIndex] * bzOld2[boundaryIndex]; // Bz at n + 1
      boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  } /* region forLoop */

  return;
}

/* Update H both in main grid and PML */
void HFieldUpdate (struct Grid *g) {
  int i,j,xStart,xStop,yStart,yStop,regionIndex,boundaryIndex;

  xStart = regionData[0].xStart;
  xStop  = regionData[0].xStop ;
  yStart = regionData[0].yStart;
  yStop  = regionData[0].yStop ;

  /* Update Hz in Main Grid */

  for (i = xStart; i < xStop; i++) {
    for (j = yStart; j < yStop; j++) {
      hz[i][j] = dahz[i][j] * hz[i][j] + dbhz[i][j] * ( ex[i][j+1] - ex[i][j] + ey[i][j] - ey[i+1][j] ); // Hz at n+1/2
    } /* jForLoop */
  } /* iForLoop */

  /* Update Hz in PML */
  boundaryIndex = 0;
  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;
    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
        hz[i][j] = hGrad1[boundaryIndex]*bz[boundaryIndex] + \
                   hGrad2[boundaryIndex]*bzOld[boundaryIndex] - \
                   hGrad3[boundaryIndex]*hz[i][j]; // Hz at n + 1/2
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  } /* region forLoop */
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
