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
#include <complex.h>
#include "empty_refl_data.h"
#include "empty_tran_data.h"
#include "fdtd_macro.h"
#include "fdtd_proto.h"


/* Implementations of Thomas Algorithm for solving Tridiagonal System
 *
 * We use this here to do the left-hand side of each field update equation for
 * a given field component and derivative direction.
 */

/* As the equations are not symmetric (they are cyclic), each update gets its
 * own function.
 * Inputs:
 * g = Grid struct that holds constants
 * i,j = int for where in the grid we are computing (either one row or one column)
 * d = Right-hand side of the update equation (d in the Thomas algorithm),
 *     size of xSize or ySize (depending on which one is needed)
 */

void exTriDiagonalSolve(struct Grid *g, int i, double *d /* d[ySize] */) {
  if (i == 0) { // Is PEC, so the field here never changes anyway
    return;
  } else {
    int k; // Thomas algorithm tracker
    double *cPrime = AllocateMemory1D(ySize, 0.0);
    double *dPrime = AllocateMemory1D(ySize, 0.0);
    // Set j = 0 term:
    cPrime[0] = c[i][0]/b[i][0];
    dPrime[0] = d[0]/b[i][0];
    for (k = 1; k < ySize; k++) {
      cPrime[k] = c[i][k]/(b[i][k] - a[i][k]*cPrime[k-1]);
      dPrime[k] = (d[k] - a[i][k]*dPrime[k-1]) / (b[i][k] - a[i][k]*cPrime[k-1]);
    }

    // Set j = ySize - 1 term:
    ex[i][ySize-1] = dPrime[ySize-1];
    for (k = ySize-2; k > -1; k--){
      ex[i][k] = dPrime[k] - cPrime[k]*ex[i][k+1];
    }

    free(cPrime);
    free(dPrime);

    return;
  }
}

void hzTriDiagonalSolve(struct Grid *g, int j, double *d /* d[xSize] */) {
  if (j == -1) { // Is PEC, so the field here never changes anyway
    return;
  } else {
    int k; // Thomas algorithm tracker
    double *cPrime = AllocateMemory1D(xSize, 0.0);
    double *dPrime = AllocateMemory1D(xSize, 0.0);
    // Set j = 0 term:
    cPrime[0] = c[0][j]/b[0][j];
    dPrime[0] = d[0]/b[0][j];
    for (k = 1; k < xSize; k++) {
      cPrime[k] = c[k][j]/(b[k][j] - a[k][j]*cPrime[k-1]);
      dPrime[k] = (d[k] - a[k][j]*dPrime[k-1]) / (b[k][j] - a[k][j]*cPrime[k-1]);
    }

    // Set i = ySize - 1 term:
    hz[xSize-1][j] = dPrime[xSize-1];
    for (k = ySize-2; k > -1; k--){
      hz[k][j] = dPrime[k] - cPrime[k]*hz[k+1][j];
    }

    free(cPrime);
    free(dPrime);

    return;
  }
}

void EFieldUpdate (struct Grid *g) {
  int i,j;
  double *d = AllocateMemory1D(ySize, 0.0);
  // See Prokopidis and Zografopoulos eq. 30


  for (i = 0; i < xSize; i++) {
    for (j = 1; j < ySize; j++) {
      exOld[i][j] = ex[i][j]; // Store previous field

      d[j] = iConst2[i][j]*ex[i][j] - ABConst[i][j]*(ex[i][j+1] - 2.0*ex[i][j] + ex[i][j-1]) + \
             ehConst[i][j]*(hz[i][j]-hz[i][j-1]) - eqConst[i][j]*qxSum[i][j];
      if(i == xSource && j == ySize - 1){
        //printf("d[254]: %.17g\n", d[254]);
      }
    }
    exTriDiagonalSolve(g,i,d);
  }

  for (i = 1; i < xSize; i++) {
    for (j = 0; j < ySize; j++) {
      eyOld[i][j] = ey[i][j]; // Store previous field

      // As Ey update depends on d/dz (which is 0 in 2D), we don't need to tri-Diagonal Solve
      // The AB term has z-derivatives which are 0 in our 2D scheme
      ey[i][j] = (iConst2[i][j]*ey[i][j] + \
                 ehConst[i][j]*(hz[i-1][j]-hz[i][j]) - \
                 eqConst[i][j]*qySum[i][j]) / (iConst1[i][j]);

                 // Might need divide by iConst1 + ABConst depending on if
                 // ey is 0 in +/- z direction or is the same as in our plane

      //eyTriDiagonalSolve(g,i,j,temp);

    } /* jForLoop */
  } /* iForLoop */

  free(d);
  return;
}

// Update Qx/Qy:
void QFieldUpdate (struct Grid *g) { // I know, it's not actually a field
  int i,j,p;

  for (p = 0; p < number_poles; p++) {
    for (i = 0; i < xSize; i++) {
      for (j = 0; j < ySize; j++) { // i=0 = pec, so don't evaluate
        // Reset sums so they don't just continuously increment:
        if (p == 0) {
          qxSum[i][j] = 0.0 + 0.0*I;
          qySum[i][j] = 0.0 + 0.0*I;
        }

        qx[p][i][j] = qConst1[p][i][j]*qx[p][i][j] + qConst2[p][i][j]*(ex[i][j]+exOld[i][j]);
        qy[p][i][j] = qConst1[p][i][j]*qy[p][i][j] + qConst2[p][i][j]*(ey[i][j]+eyOld[i][j]);
        // Accumulate sums
        qxSum[i][j] += creal(qSumC[p][i][j]*qx[p][i][j]);
        qySum[i][j] += creal(qSumC[p][i][j]*qy[p][i][j]);
      } /* jForLoop */
    } /* iForLoop */
  } /* pForLoop */

  return;
}


/* Update H both in main grid and PML */
void HFieldUpdate (struct Grid *g) {
  int i,j;
  double *d = AllocateMemory1D(xSize, 0.0);

  for (j = 0; j < ySize; j++) { // THE ORDER HERE IS BACKWARDS. There might be a better way as this is an inefficient way to access these elements
    for (i = 1; i < xSize; i++) {
      d[i] = hz[i][j] - ABConst[i][j]*(hz[i+1][j] - 2*hz[i][j] + hz[i-1][j]) + \
             heConst[i][j]*((ex[i][j+1] - ex[i][j]) - (ey[i+1][j] - ey[i][j])); // INCORRECT YEE LATTICE HERE MAYBE (maybe not...?)
             // This also assume magnetic permeability = 1
    } /* iForLoop */
    hzTriDiagonalSolve(g,j,d); // Now do the tri-diagonal solving
  } /* jForLoop */

  free(d);
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
