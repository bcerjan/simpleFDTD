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
#include "fdtd_macro.h"
#include "fdtd_proto.h"
#include "structure_funcs.h"
#include "ezinc.h"
#include "material_data.h"
#include "array_proto.h"



/* Functions to get material data (or create it) from our stored file */
double getMatPermittivity(int metalChoice, double objectIndex) {
  if( metalChoice > -1 ) {
    return materialData[metalChoice].epsInf;
  } else {
    return objectIndex*objectIndex;
  }
}
double getMatConductivity(int metalChoice) {
  if( metalChoice > -1 ) {
    return materialData[metalChoice].conductivity;
  } else {
    return 0.0;
  }
}
double getMatPermeability(int metalChoice) {
  if( metalChoice > -1 ) {
    return materialData[metalChoice].permeability;
  } else {
    return 1.0;
  }
}
double getMatResistivity(int metalChoice) {
  if( metalChoice > -1 ) {
    return 1.0/materialData[metalChoice].conductivity;
  } else {
    return 0.0;
  }
}

// Function to get the index we're using to reference our "empty" runs by
// this is, unfortunately, the index of the (refractive) index.
int getIndexIndex(double environmentIndex) {

  return lround(environmentIndex*10.0) - 10;
}
/* Function to initialize our Grid and it's associated constants */

/* Inputs:
struct Grid *g : The grid struct we are working with currently
int metalChoice : 0 = Al, 1 = Au, 2 = Ag, 3 = Cu, other = SiO2
int objectChoice : 0 = Disk, 1 = Rectangle, 2 = triangle, other = no structure
double objectSize : Object size (in nm)
double environmentIndex : Refractive Index of the environment
*/

void  InitializeFdtd (struct Grid *g, int metalChoice, int objectChoice,
  double objectXSize, double objectYSize, double environmentIndex, double objectIndex )
{

    double  mediaPermittivity[MEDIACONSTANT] = {environmentIndex*environmentIndex, getMatPermittivity(metalChoice, objectIndex)};    // eps, index=0 is for vacuum, index=1 is for the metallic cylinder
    double  mediaConductivity[MEDIACONSTANT] = {0.0, getMatConductivity(metalChoice)}; // sig,
    double  mediaPermeability[MEDIACONSTANT] = {1.0, getMatPermeability(metalChoice)};    // mur
    double  mediaResistivity[MEDIACONSTANT] = {0.0, getMatResistivity(metalChoice)};     // sim
    complex double  mediaA[MEDIACONSTANT][MAX_POLES] = {{0.0}};
    complex double  mediaC[MEDIACONSTANT][MAX_POLES] = {{0.0}};
    double  magneticPermeability0,electricalPermittivity0,frequency,wavelength,angularFrequency;
    double  temporary,electricalImpedance0;
    int  i,j,p,media,xSizeMain,ySizeMain,numFreqs;
    int  xCenter,yCenter;

    /***********************************************************************/
    //     Fundamental constants
    /***********************************************************************/

    pi  = (acos(-1.0));
    speedOfLight = 2.99792458e8;                          //speed of light in free space (meters/second)
    magneticPermeability0 = 4.0 * pi * 1.0e-7;            //permeability of free space
    electricalPermittivity0 = 1.0 / (speedOfLight * speedOfLight * magneticPermeability0);       //permittivity of free space
    electricalImpedance0 = 377.0;                         // Electrical impedance of free space (377 Ohm) (eta-0)

    // 5.0e14 is ~600 nm wavelength
    frequency = 5.0e+14;                                  //center frequency of source excitation (Hz)
    wavelength = speedOfLight / frequency  ;              //center wavelength of source excitation
    angularFrequency = 2.0 * pi * frequency;              //center frequency in radians

    // Set grid spacing:
    dx = 10.0e-9; // 10 nm
    double dxnm = dx*1e9; // Grid step size in nm

    //courantS = 1.0/2.0;
    courantS = 2.5/2.0;
    dt = courantS * dx / speedOfLight;

    // Divide by environment index. This means that if the object goes through the boundary, weird stuff might happen
    absConst = ( (speedOfLight/environmentIndex)*dt - dx) / ((speedOfLight/environmentIndex)*dt + dx);
//printf( "dx: %f\n", dx );
//printf( "dt: %f\n", dt );
//printf( "courantS: %f\n", courantS );
    /***********************************************************************/
    //     Grid parameters
    /***********************************************************************/

    xSizeMain = 300; //300                             // number of main grid cells in x-direction
    ySizeMain = 250; //250                              // number of main grid cells in y-direction
    //abcSize = ABCSIZECONSTANT;                    // thickness of PML region
    xSize = xSizeMain; //+ 2 * abcSize;              // number of total grid cells in x-direction
    ySize = ySizeMain; //+ 2 * abcSize;              // number of total grid cells in y-direction

    xSource = 20;                          //location of z-directed hard source

    envIndex = environmentIndex;                  // Background refractive index

    maximumIteration = NUMBEROFITERATIONCONSTANT;                 //total number of time steps



    /***********************************************************************/
    //     Material parameters
    /***********************************************************************/

    media = MEDIACONSTANT;        // number of different medias, ie 2: vacuum, metallicCylinder

    refractiveIndexIndex = getIndexIndex(environmentIndex);
    //printf("refractiveIndexIndex %i\n", refractiveIndexIndex);
    // Number of poles in our dielectric function
    if (metalChoice > -1 && objectChoice > -1) {
      number_poles = materialData[metalChoice].num_poles;
      for (p = 0; p < number_poles; p++) {
        mediaA[1][p] = materialData[metalChoice].params[p].ap;
        mediaC[1][p] = materialData[metalChoice].params[p].cp;
      }
    } else {
      number_poles = 0;
    }


    /***********************************************************************/
    //     Wave excitation
    /***********************************************************************/

    ezIncInit(g, environmentIndex);

    /***********************************************************************/
    //     Field arrays
    /***********************************************************************/

    ex = AllocateMemory(xSize,    ySize + 1, 0.0 );        // 1 extra in y direction for pec
    ey = AllocateMemory(xSize + 1,ySize,     0.0 );        // 1 extra in x direction for pec
    hz = AllocateMemory(xSize + 1,ySize + 1, 0.0 );

    e2Field = AllocateMemory(xSize, ySize, 0.0);           // E^2 for plotting
    edgeMat = AllocateMemory(xSize, ySize, 0.0);

    /* Polarization Current Fields */
    complex double initArray[MAX_POLES] = {0.0};
    qx = AllocateComplexMemory3D(number_poles, xSize, ySize, initArray);
    qy = AllocateComplexMemory3D(number_poles, xSize, ySize, initArray);

    exOld = AllocateMemory(xSize, ySize + 1, 0.0);
    eyOld = AllocateMemory(xSize + 1, ySize, 0.0);
    hzOld = AllocateMemory(xSize + 1, ySize + 1, 0.0);


    /*printf("cjjTemp[0]: %.5e\n", cjjTemp[0]);
    printf("cjjTemp[1]: %.5e\n", cjjTemp[1]);
    printf("cjeTemp[0]: %.5e\n", cjeTemp[0]);
    printf("cjeTemp[1]: %.5e\n", cjeTemp[1]);*/

    /***********************************************************************/
    //     DFT Array Initialization
    /***********************************************************************/
    numFreqs = NUMBERDFTFREQS;
    reflDFT = AllocateMemory1D(numFreqs, 0.0);
    tranDFT = AllocateMemory1D(numFreqs, 0.0);

    reEyReflDFT = AllocateMemory(numFreqs, ySize, 0.0);
    imEyReflDFT = AllocateMemory(numFreqs, ySize, 0.0);
    reEyTranDFT = AllocateMemory(numFreqs, ySize, 0.0);
    imEyTranDFT = AllocateMemory(numFreqs, ySize, 0.0);

    reHzReflDFT = AllocateMemory(numFreqs, ySize, 0.0);
    imHzReflDFT = AllocateMemory(numFreqs, ySize, 0.0);
    reHzTranDFT = AllocateMemory(numFreqs, ySize, 0.0);
    imHzTranDFT = AllocateMemory(numFreqs, ySize, 0.0);

    reflXPos = 5;
    tranXPos = xSize - 5;

    double waveMin, waveMax; // Min / Max wavelength in PPW
    waveMin = (400.0e-9) / dx; // 400 nm
    waveMax = (800.0e-9) / dx; // to 800 nm

    double maxDFTTime = DFTPADDEDTIME; // Since we 0-pad the DFT we have a ficticious time that is longer than the simulation
    // Find frequencies:
    for (i = 0; i < numFreqs; i++) {
      // Calculate evenly spaced frequencies from 400 to 800 nm
      temporary = waveMin + ((double  )i * (waveMax - waveMin) / ((double  )numFreqs - 1));
      kList[i] = (int  )round(maxDFTTime * courantS / temporary); // Schneider 5.29

      // Now convert back to find the "real" wavelength we are working with:
      wavelengthList[i] = maxDFTTime * courantS / (double  )kList[i]; // Same, but inverted for wavelength
      //printf("kList[%i]: %i\n", i, kList[i] );
      //printf("wavelengthList[%i]: %f\n", i, wavelengthList[i] );
    } /* iForLoop */

    /***********************************************************************/
    //     Media coefficients
    /***********************************************************************/

    heConst = AllocateMemory(xSize, ySize, dt/(dx*magneticPermeability0*mediaPermeability[0]));
    ehConst = AllocateMemory(xSize, ySize, dt/(dx*electricalPermittivity0*mediaPermittivity[0]));
    eqConst = AllocateMemory(xSize, ySize, 4.0*dt/(electricalPermittivity0*mediaPermittivity[0]));
    qxSum   = AllocateMemory(xSize, ySize, 0.0);
    qySum   = AllocateMemory(xSize, ySize, 0.0);
    qConst1 = AllocateComplexMemory3D(number_poles, xSize, ySize, initArray);
    qConst2 = AllocateComplexMemory3D(number_poles, xSize, ySize, initArray);
    qSumC   = AllocateComplexMemory3D(number_poles, xSize, ySize, initArray);
    ABConst = AllocateMemory(xSize, ySize, dt*dt/(4.0*magneticPermeability0*electricalPermittivity0*dx*dx*mediaPermeability[0]*mediaPermittivity[0]));


    /*printf("heConst: %.5e\n", heConst[0][0]);
    printf("ehConst: %.5e\n", ehConst[0][0]);
    printf("eqConst: %.5e\n", eqConst[0][0]);
    printf("ABConst: %.5e\n", ABConst[0][0]);*/

    // Temp storage values:
    complex double **qC1 = AllocateComplexMemory(media, number_poles, 0.0);
    complex double **qC2 = AllocateComplexMemory(media, number_poles, 0.0);
    double iCSum[MEDIACONSTANT] = {0.0};
    double iC1[MEDIACONSTANT] = {0.0};
    double iC2[MEDIACONSTANT] = {0.0};
    /* For terms and formulation, see: A Unified FDTD/PML Scheme Based on Critical
      Points for Accurate Studies of Plasmonic Structures. Prokopidis and
      Zografopoulos, JOURNAL OF LIGHTWAVE TECHNOLOGY, VOL. 31, NO. 15, AUGUST 1, 2013 */

    // Electric / Polarization field update constants:
    // Critical Points terms:
    for (i = 0; i < media; i++) {
      for (p = 0; p < number_poles; p++){
        /* C-P Parameters */
        qC1[i][p] = (2.0 + mediaA[i][p]*dt) / (2.0 - mediaA[i][p]*dt);
        qC2[i][p] = electricalPermittivity0*mediaC[i][p]*dt / (2.0 - mediaA[i][p]*dt);

        iCSum[i] += creal( mediaC[i][p]/(2.0 - mediaA[i][p]*dt) );
      } /* pForLoop */
    } /* iForLoop */

    // Set Q factors to vacuum values everywhere:
    for (p = 0; p < number_poles; p++) {
      for (i = 0; i < xSize; i++) {
        for (j = 0; j < ySize; j++) {
          qConst1[p][i][j] = qC1[0][p];
          qConst2[p][i][j] = qC2[0][p];

          qSumC[p][i][j] = mediaA[0][p] / (2.0-mediaA[0][p]*dt);
        } /* jForLoop */
      } /* iForLoop */
    } /* pForLoop */

    // Identity matrix coefficients
    for (i = 0; i < media; i++) {
      iC1[i] = 1.0 + \
               mediaConductivity[i]*dt/(2.0*electricalPermittivity0*mediaPermittivity[i]) + \
               2.0*dt*iCSum[i]/(mediaPermittivity[i]);
      iC2[i] = 1.0 - \
               mediaConductivity[i]*dt/(2.0*electricalPermittivity0*mediaPermittivity[i]) - \
               2.0*dt*iCSum[i]/(mediaPermittivity[i]);
    } /* iForLoop */

    iConst1 = AllocateMemory(xSize, ySize, iC1[0]);
    iConst2 = AllocateMemory(xSize, ySize, iC2[0]);

    /*printf("iC1[0]: %.5e\n", iC1[0]);
    printf("iC1[1]: %.5e\n", iC1[1]);
    printf("iC2[0]: %.5e\n", iC2[0]);
    printf("iC2[1]: %.5e\n", iC2[1]);*/

    /***********************************************************************/
    //     Grid Coefficients
    /***********************************************************************/

    /* Array to track where our object is/is not */
    object_locs = AllocateMemory(xSize, ySize, 0.0); // This really shouldn't be an array of doubles -- we're wasting memory here

    // Initialize structure functions:
    xCenter = xSize - 30; // In grid units
    yCenter = ySize / 2; // ""
    structInit(xCenter, yCenter);
printf("Strucutre Init...\n" );

    // Sanity checks on input sizes, if statements make sure we're never trying
    // to update outside the grid or beyond the monitor positions:
    double x_size,y_size,x_steps,y_steps;
    if( objectXSize < 0.0 ) {
      x_size = 0.0;
      objectXMax = 2; // Should never be used, but just in case
      objectXMin = 1; // ""
    } else {
      x_size = objectXSize * dx / dxnm;
      x_steps = objectXSize * 1.0e-9 / dx; // Assumes units are in nm
      objectXMax = xCenter + (int )ceil(x_steps/2.0);
      if (objectXMax > 2*xSize-tranXPos+1) { objectXMax = 2*xSize-tranXPos+1; }
      objectXMin = xCenter - (int )floor(x_steps/2.0);
      if (objectXMin < reflXPos+1) { objectYMax = reflXPos+1; }
    } /* ifBlock */

    if( objectYSize < 0.0 ) {
      y_size = 0.0;
      objectYMax = 2; // Should never be used, but just in case
      objectYMin = 1; // ""
    } else {
      y_size = objectYSize * dx / dxnm;
      y_steps = objectYSize * 1.0e-9 / dx; // Assumes units are in nm
      objectYMax = yCenter + (int )ceil(y_steps/2.0);
      if (objectYMax > ySize-1) { objectYMax = ySize-1; }
      objectYMin = yCenter - (int )floor(y_steps/2.0);
      if (objectYMin < 1) { objectYMin = 1; }
    } /* ifBlock */

    // Switch Block to pick structure geometry (default is no object):
    switch (objectChoice) {
      case 0: // Disk
        addDisk(g, x_size/2.0, y_size/2.0); // divide by 2 to get radius
        break;

      case 1: // Block
        addRect(g, x_size, y_size);
        break;

      case 2: // Triangle
        addTriangle(g, x_size, y_size);
        break;
    } /* switch */

    // Add structure (assumes only one material besides background):
    for (i = 0; i < xSize; i++) {
      for (j = 0; j < ySize; j++) {
        if (object_locs[i][j] > 0.5) {
          // Material Constants:
          ABConst[i][j] = dt*dt/(4.0*magneticPermeability0*electricalPermittivity0*dx*dx*mediaPermeability[1]*mediaPermittivity[1]);
          ehConst[i][j] = dt/(dx*electricalPermittivity0*mediaPermittivity[1]);
          eqConst[i][j] = 4.0*dt/(electricalPermittivity0*mediaPermittivity[1]);
          heConst[i][j] = dt/(dx*magneticPermeability0*mediaPermeability[1]);

          iConst1[i][j] = iC1[1];
          iConst2[i][j] = iC2[1];

           /* C-P C coefficients */
          for(p = 0; p < number_poles; p++) {
            /* Polarization Constants: */
            qConst1[p][i][j] = qC1[1][p];
            qConst2[p][i][j] = qC2[1][p];

            qSumC[p][i][j] = mediaA[1][p] / (2.0-mediaA[1][p]*dt);
          } /* pForLoop */
        } /* ifBlock */
      } /* jForLoop */
    } /* iForLoop */

    freeComplexDoublePtr(qC1, media);
    freeComplexDoublePtr(qC2, media);

printf("Structure Added...\n" );
printf("Structure Choice: %i\n", objectChoice);
printf("objectXSize: %f\n", objectXSize);
printf("object_locs[xCent][yCent]: %f\n", object_locs[xCenter][yCenter]);

    /***********************************************************************/
    //     Tridiagonal Solver Coefficients
    /***********************************************************************/
    // This needs to be done last as they take the values from above

    // For dx^2 terms (used in hzTriDiagonalSolve):
    ahz = AllocateMemory(xSize, ySize, 0.0);
    bhz = AllocateMemory(xSize, ySize, 0.0);
    chz = AllocateMemory(xSize, ySize, 0.0);
    // for dy^2 terms (used in exTriDiagonalSolve):
    aex = AllocateMemory(xSize, ySize, 0.0);
    bex = AllocateMemory(xSize, ySize, 0.0);
    cex = AllocateMemory(xSize, ySize, 0.0);
    // No dz^2 terms as this is 2D!

    for (i = 0; i < xSize; i++) {
      for (j = 0; j < ySize; j++) {
        ahz[i][j] = -1.0 * ABConst[i][j];
        //bhz[i][j] = iConst1[i][j] + 2.0*ABConst[i][j];
        bhz[i][j] = 1.0 + 2.0*ABConst[i][j]; // no e(w) term for Hz update
        chz[i][j] = ahz[i][j];

        aex[i][j] = -1.0 * ABConst[i][j];
        bex[i][j] = iConst1[i][j] + 2.0*ABConst[i][j];
        cex[i][j] = aex[i][j];
      }
    }

    /***********************************************************************/
    //     Mur ABC Setup
    /***********************************************************************/



    // Add Mur ABC conditions from "Mur Absorbing Boundary Condition for 2-D
    // Leapgfrog ADI-FDTD Method", Gan and Tan, IEEE Asia-Pacific Conf. 2012

    // As Ex depends on j-terms:
    /* Here the Ex ABC is implemented at j = 1 instead of j = 0. This is because
       (emprically) it produces instabilities of it's at 0. I believe the issue
       stems from the staggering of the nodes (i.e. the Yee cell) and so if you
       put the ABC at j = 0, Ex grows non-physically at j = 1 (as Hz -> 0 at
       j = 0, but not at j = 1). The Update equation for Ex also ignores the
       j = 0 row. */

    for (i = 0; i < xSize; i++) {
      aex[i][1] = 0.0;
      bex[i][1] = 1.0;
      cex[i][1] = -1.0*absConst;

      aex[i][ySize-1] = -1.0*absConst;
      bex[i][ySize-1] = 1.0;
      cex[i][ySize-1] = 0.0;
    }

    // While Hz depends on i-terms:
    for (j = 0; j < ySize; j++) {
      ahz[0][j] = 0.0;
      bhz[0][j] = 1.0;
      chz[0][j] = -1.0*absConst;

      ahz[xSize-1][j] = -1.0*absConst;
      bhz[xSize-1][j] = 1.0;
      chz[xSize-1][j] = 0.0;
    }

/*    printf("ahz: %.17g\n", ahz[4][0]);
    printf("bhz: %.17g\n", bhz[4][0]);

    printf("aex: %.17g\n", aex[4][0]);
    printf("bex: %.17g\n", bex[4][0]);
*/

    // all done with Initialization!

    return;

}
