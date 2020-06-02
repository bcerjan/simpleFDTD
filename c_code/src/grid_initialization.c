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
    double  drudePlasma[MEDIACONSTANT] = {0.0};
    double  drudeDamping[MEDIACONSTANT] = {HUGE_VAL}; // initialize vacuum to have huge damping
    double  mediaCa[MEDIACONSTANT];
    double  mediaCb[MEDIACONSTANT];
    double  mediaDa[MEDIACONSTANT];
    double  mediaDb[MEDIACONSTANT];
    double  magneticPermeability0,electricalPermittivity0,frequency,wavelength,angularFrequency;
    double  reflectionCoefficient0,gradingOrder,temporary,electricalImpedance0,temp1,temp2,temp3;
    int  i,j,k,p, boundaryDataSize, media, boundaryIndex,xSizeMain,ySizeMain,numFreqs;
    int  abcSize ;
    double  temporaryi,temporaryj,distance2 ;
    int  xCenter,yCenter;
    double  x,x1,x2;
    double  electricalConductivityMaximum, kappaMaximum, boundaryWidth, alphaPML;
    double  gradientConductivity, gradientResistivity, boundaryFactor, denominator;
    double  gradientCa1[ABCSIZECONSTANT];
    double  gradientCb1[ABCSIZECONSTANT];
    double  gradientDa1[ABCSIZECONSTANT];
    double  gradientDb1[ABCSIZECONSTANT];

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
    dt = dx / (2.0 * speedOfLight);
    courantS = (dt * speedOfLight) / (dx); // Changed from original -- might break things?
//printf( "dx: %f\n", dx );
//printf( "dt: %f\n", dt );
//printf( "courantS: %f\n", courantS );
    /***********************************************************************/
    //     Grid parameters
    /***********************************************************************/

    xSizeMain = 300;                              // number of main grid cells in x-direction
    ySizeMain = 250;                               // number of main grid cells in y-direction
    abcSize = ABCSIZECONSTANT;                    // thickness of PML region
    xSize = xSizeMain + 2 * abcSize;              // number of total grid cells in x-direction
    ySize = ySizeMain + 2 * abcSize;              // number of total grid cells in y-direction

    boundaryDataSize  = 2 * xSize * abcSize;                      // front edge + back edge
    boundaryDataSize += 2 * (abcSize * (ySize - 2 * abcSize));    // left + right edges

    //xSource = 50 + abcSize;                          //location of z-directed hard source
    //ySource = 50 + abcSize;                          //location of z-directed hard source
    xSource = 20 + abcSize;                       //location of z-directed hard source
    ySource = ySize / 2;                          //location of z-directed hard source

    envIndex = environmentIndex;                  // Background refractive index

    maximumIteration = NUMBEROFITERATIONCONSTANT;                 //total number of time steps

    reflectionCoefficient0 = 1.0e-7;              // for PML, Nikolova part4 p.25
    gradingOrder = 3;                             // for PML, (m) was 2;  optimal values: 2 <= m <= 6,  Nikolova part4 p.29

    /***********************************************************************/
    //     Material parameters
    /***********************************************************************/

    media = MEDIACONSTANT;        // number of different medias, ie 2: vacuum, metallicCylinder

    refractiveIndexIndex = getIndexIndex(environmentIndex);
    printf("refractiveIndexIndex %i\n", refractiveIndexIndex);
    // Number of poles in our dielectric function
    if (metalChoice > -1 && objectChoice > -1) {
      number_poles = materialData[metalChoice].num_poles;
    } else {
      number_poles = 0;
    }

    // Constants for our approximation:
    double **Ap,**Omegap,**phip,**Gammap;
    Ap = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);
    Omegap = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);
    phip = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);
    Gammap = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);

    for (i = 1; i < media; i++) {
      for (p = 0; p < number_poles; p++) {
        Ap[i][p] = materialData[metalChoice].params[p].bigA;
        Omegap[i][p] = materialData[metalChoice].params[p].Omega;
        phip[i][p] = materialData[metalChoice].params[p].phi;
        Gammap[i][p] = materialData[metalChoice].params[p].Gamma;
      } /* pForLoop */
    } /* iForLoop */

    for (i = 1; i < media; i++) {
      drudePlasma[i] = materialData[metalChoice].drudePlasma;
      drudeDamping[i] = materialData[metalChoice].drudeDamping;
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
    hz = AllocateMemory(xSize,    ySize,     0.0 );
    hzy = AllocateMemory1D(boundaryDataSize, 0.0);         // For split-field data in PML

    hzOld = AllocateMemory1D(boundaryDataSize, 0.0);

    e2Field = AllocateMemory(xSize, ySize, 0.0);           // E^2 for plotting
    edgeMat = AllocateMemory(xSize, ySize, 0.0);

    /* Polarization Current Fields */
    double initArray[MAX_POLES] = {0.0};
    px = AllocateMemory3D(number_poles, xSize, ySize, initArray);
    py = AllocateMemory3D(number_poles, xSize, ySize, initArray);
    pxOld = AllocateMemory3D(number_poles, xSize, ySize, initArray);
    pyOld = AllocateMemory3D(number_poles, xSize, ySize, initArray);
    pxOld2 = AllocateMemory3D(number_poles, xSize, ySize, initArray);
    pyOld2 = AllocateMemory3D(number_poles, xSize, ySize, initArray);
    exOld = AllocateMemory(xSize, ySize + 1, 0.0);
    exOld2 = AllocateMemory(xSize, ySize + 1, 0.0);
    eyOld = AllocateMemory(xSize + 1, ySize, 0.0);
    eyOld2 = AllocateMemory(xSize + 1, ySize, 0.0);

    pxDrude = AllocateMemory(xSize, ySize, 0.0);
    pyDrude = AllocateMemory(xSize, ySize, 0.0);

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

    reflXPos = abcSize + 5;
    tranXPos = xSize - abcSize - 5;

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

    double **c1,**c2,**c3,**c4,**c5,*c3TempSum,*c4TempSum,*c5TempSum;
    c1 = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);
    c2 = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);
    c3 = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);
    c4 = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);
    c5 = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);
    c3TempSum = AllocateMemory1D(MEDIACONSTANT, 0.0);
    c4TempSum = AllocateMemory1D(MEDIACONSTANT, 0.0);
    c5TempSum = AllocateMemory1D(MEDIACONSTANT, 0.0);


    /* Pre-compute values that will be placed in to cjj and cje */
    double cjjTemp[media], cjeTemp[media];
    double AZero,AOne,BZero,BOne,BTwo,cP,d2,cond;
    double nDamping;


    /* For terms and formulation, see: A Unified FDTD/PML Scheme Based on Critical
      Points for Accurate Studies of Plasmonic Structures. Prokopidis and
      Zografopoulos, JOURNAL OF LIGHTWAVE TECHNOLOGY, VOL. 31, NO. 15, AUGUST 1, 2013 */

    // Electric / Polarization field update constants:
    // Critical Points terms:
    for (i = 0; i < media; i++) {
      for (p = 0; p < number_poles; p++){
        /* C-P Parameters */
        AZero = 2.0*electricalPermittivity0*Ap[i][p]*Omegap[i][p] * \
          ( Omegap[i][p]*cos(phip[i][p]) - Gammap[i][p]*sin(phip[i][p]) );
        AOne = -2.0*electricalPermittivity0*Ap[i][p]*Omegap[i][p]*sin(phip[i][p]);
        BZero = pow(Omegap[i][p],2) + pow(Gammap[i][p],2);
        BOne = 2.0*Gammap[i][p];
        BTwo = 1.0;
        cP = BTwo/(dt*dt) + BOne/(2.0*dt) + BZero/4.0;

        c1[i][p] = (2.0*BTwo/(dt*dt) - BZero/2.0)/cP;
        c2[i][p] = (BOne/(2.0*dt) - BTwo/(dt*dt) - BZero/4.0)/cP;
        c3[i][p] = (AZero/4.0 + AOne/(2.0*dt))/cP;
        c4[i][p] = AZero/(2.0*cP);
        c5[i][p] = (AZero/4.0 - AOne/(2.0*dt))/cP;

        c3TempSum[i] += c3[i][p];
        c4TempSum[i] -= c4[i][p];
        c5TempSum[i] += c5[i][p];
      } /* pForLoop */

      /* Drude d2 parameter: */
      d2 = electricalPermittivity0*drudePlasma[i]*drudePlasma[i]*dt / \
        ( (2.0 + drudeDamping[i]*dt)*drudeDamping[i] );
      cond = electricalPermittivity0*drudePlasma[i]*drudePlasma[i] / drudeDamping[i];

      // See eq. 23 in the paper:
      c3TempSum[i] += electricalPermittivity0*mediaPermittivity[i] + \
        0.5*cond*dt - d2;
      c4TempSum[i] += electricalPermittivity0*mediaPermittivity[i] - \
        0.5*cond*dt + d2;
    } /* iForLoop */

    // Field Update Constants:
    for (i = 0; i < media; i++) {
      //mediaCa[i] = 1.0/c3TempSum[i]; // Multiplied by E and P terms in E update equation
      mediaCa[i] = 1.0; // Only used in PML
      mediaCb[i] = dt / (dx * c3TempSum[i]); // Multiplied by H in E update
      //temporary = mediaConductivity[i] * dt / (2.0 * electricalPermittivity0 * mediaPermittivity[i]); // Schneider 8.13 - 8.18
      mediaDa[i] = 1.0; // assuming magnetic conductivity is 0
      mediaDb[i] = dt / ( dx*magneticPermeability0 * mediaPermeability[i]); // Multiplied by E in H update equation
      //mediaDb[i] = dt;
    } /* iForLoop */

    /***********************************************************************/
    //     Grid Coefficients
    /***********************************************************************/
    /*printf("c3TempSum[0]: %.17g\n", c3TempSum[0]);
    printf("c3TempSum[1]: %.17g\n", c3TempSum[1]);
    printf("c4TempSum[0]: %.17g\n", c4TempSum[0]);
    printf("c4TempSum[1]: %.17g\n", c4TempSum[1]);
    printf("c5TempSum[0]: %.17g\n", c5TempSum[0]);
    printf("c5TempSum[1]: %.17g\n", c5TempSum[1]);
    printf("mediaDa[0]: %f\n", mediaDa[0]);
    printf("mediaDa[1]: %f\n", mediaDa[1]);
    printf("mediaDb[0]: %f\n", mediaDb[0]);
    printf("mediaDb[1]: %f\n", mediaDb[1]);
    printf("mediaCa[0]: %f\n", mediaCa[0]);
    printf("mediaCa[1]: %f\n", mediaCa[1]);
    printf("mediaCb[0]: %f\n", mediaCb[0]);
    printf("mediaCb[1]: %f\n", mediaCb[1]);*/
    //     Initialize entire grid to free space

    // Polarization grid values:
    c1SumX = AllocateMemory(xSize, ySize, 0.0); // No temp sums for c1 and c2, and this is 0 anyway
    c2SumX = AllocateMemory(xSize, ySize, 0.0);
    c1SumY = AllocateMemory(xSize, ySize, 0.0);
    c2SumY = AllocateMemory(xSize, ySize, 0.0);
    c3Sum = AllocateMemory(xSize, ySize, c3TempSum[0]);
    c4Sum = AllocateMemory(xSize, ySize, c4TempSum[0]);
    c5Sum = AllocateMemory(xSize, ySize, c5TempSum[0]);


    // Make so we don't have to loop over each array twice on initialization...
    c1Grid = AllocateMemory3D(number_poles, xSize, ySize, c1[0]);
    c2Grid = AllocateMemory3D(number_poles, xSize, ySize, c2[0]);
    c3Grid = AllocateMemory3D(number_poles, xSize, ySize, c3[0]);
    c4Grid = AllocateMemory3D(number_poles, xSize, ySize, c4[0]);
    c5Grid = AllocateMemory3D(number_poles, xSize, ySize, c5[0]);

    d1Grid = AllocateMemory(xSize, ySize, 0.0);
    d2Grid = AllocateMemory(xSize, ySize, 0.0);

/*    for (i = 0; i < xSize; i++) {
      for (j = 0; j < ySize; j++) {
        for (p = 0; p < number_poles; p++) {
          c1Grid[p][i][j] = c1[0][p];
          c2Grid[p][i][j] = c2[0][p];
          c3Grid[p][i][j] = c3[0][p];
          c4Grid[p][i][j] = c4[0][p];
          c5Grid[p][i][j] = c5[0][p];
        }*/ /* pForLoop */
      //} /* jForLoop */
    //} /* iForLoop */



    /* Values for PML updates: */
    caex = AllocateMemory(xSize, ySize, mediaCa[0] );     // note: don't need to allocate for pec region, as it is not evaluated
    cbex = AllocateMemory(xSize, ySize, mediaCb[0] );     // also: Initialize the entire grid to vacuum.
    caey = AllocateMemory(xSize, ySize, mediaCa[0] );
    cbey = AllocateMemory(xSize, ySize, mediaCb[0] );
    dahz = AllocateMemory(xSize, ySize, mediaDa[0] );
    dbhz = AllocateMemory(xSize, ySize, mediaDb[0] );
    dahzy = AllocateMemory1D(boundaryDataSize, mediaDa[0] );        // for the split-field data for hz in the pml regions
    dbhzy = AllocateMemory1D(boundaryDataSize, mediaDb[0] );        // for the split-field data for hz in the pml regions

    /* Array to track where our object is/is not */
    object_locs = AllocateMemory(xSize, ySize, 0.0); // This really shouldn't be an array of doubles -- we're wasting memory here

    // Initialize structure functions:
    xCenter = xSize - abcSize - 30; // In grid units
    yCenter = ySize / 2; // ""
    structInit(xCenter, yCenter);
printf("Strucutre Init...\n" );

    // Sanity checks on input sizes:
    double x_size,y_size;
    if( objectXSize < 0.0 ) {
      x_size = 0.0;
    } else {
      x_size = objectXSize * dx / dxnm;
    } /* ifBlock */

    if( objectYSize < 0.0 ) {
      y_size = 0.0;
    } else {
      y_size = objectYSize * dx / dxnm;
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
          c3Sum[i][j] = c3TempSum[1];
          c4Sum[i][j] = c4TempSum[1];
          c5Sum[i][j] = c5TempSum[1];

          caex[i][j] = mediaCa[1];
          caey[i][j] = mediaCa[1];
          cbex[i][j] = mediaCb[1];
          cbey[i][j] = mediaCb[1];
          dahz[i][j] = mediaDa[1];
          dbhz[i][j] = mediaDb[1];
           /* C-P C coefficients */
          for(p = 0; p < number_poles; p++) {
            /* Polarization Constants: */
            c1Grid[p][i][j] = c1[1][p];
            c2Grid[p][i][j] = c2[1][p];
            c3Grid[p][i][j] = c3[1][p];
            c4Grid[p][i][j] = c4[1][p];
            c5Grid[p][i][j] = c5[1][p];
          } /* pForLoop */

          /* Drude Coefficients */
          d1Grid[i][j] = (2.0 - drudeDamping[1]*dt)/(2.0 + drudeDamping[1]*dt);
          d2Grid[i][j] = electricalPermittivity0*drudePlasma[1]*drudePlasma[1]*dt / \
            ( (2.0 + drudeDamping[1]*dt)*drudeDamping[1] );
        } /* ifBlock */
      } /* jForLoop */
    } /* iForLoop */

    // Free arrays only used here:
    freeDoublePtr(c1, MEDIACONSTANT);
    freeDoublePtr(c2, MEDIACONSTANT);
    freeDoublePtr(c3, MEDIACONSTANT);
    freeDoublePtr(c4, MEDIACONSTANT);
    freeDoublePtr(c5, MEDIACONSTANT);
    freeDoublePtr(Ap, MEDIACONSTANT);
    freeDoublePtr(Gammap, MEDIACONSTANT);
    freeDoublePtr(phip, MEDIACONSTANT);
    freeDoublePtr(Omegap, MEDIACONSTANT);
    free(c3TempSum);
    free(c4TempSum);
    free(c5TempSum);

printf("Strucutre Added...\n" );
    /***********************************************************************/
    //     Initialize the RegionDataValues structure
    /***********************************************************************/

    // regions are arranged in this order: center(main), front, back, left, right
    // note: the region order is "important" in order to make the split-field data line up right for hzy (see the Fill the PML section below)
    // this structure is for calculating Hz. (need to break Hz up into pml (split-field) regions and main grid)
    // how they are used: loopIndex = regionData.start, loopIndex < regionData.stop
    regionData[0].xStart = abcSize;                    // main grid
    regionData[0].xStop  = abcSize + xSizeMain;
    regionData[0].yStart = abcSize;
    regionData[0].yStop  = abcSize + ySizeMain;
    regionData[1].xStart = 0;                          // bottom grid
    regionData[1].xStop  = xSize;
    regionData[1].yStart = 0;
    regionData[1].yStop  = abcSize;
    regionData[2].xStart = 0;                          // top grid
    regionData[2].xStop  = xSize;
    regionData[2].yStart = ySize - abcSize;
    regionData[2].yStop  = ySize;
    regionData[3].xStart = 0;                          // left grid
    regionData[3].xStop  = abcSize;
    regionData[3].yStart = abcSize;
    regionData[3].yStop  = abcSize + ySizeMain;
    regionData[4].xStart = xSize - abcSize;            // right grid
    regionData[4].xStop  = xSize;
    regionData[4].yStart = abcSize;
    regionData[4].yStop  = abcSize + ySizeMain;

    /***********************************************************************/
    //     Fill the PML regions    ---  (Caution...Here there be Tygers!)
    /***********************************************************************/

    // The most important part of the PML fdtd simulation is getting the
    // PML Coefficients correct. Which requires getting the correct PML gradient and
    // positioning the coefficients correctly on the x-y grid.

    // ALERT: It is possible to make a mistake here, and yet the simulation may appear to
    // be working properly. However a detailed analysis of reflections off the PML
    // will show they may be (much) larger than those for a correctly designed PML.

    boundaryWidth = (double  )abcSize * dx;    // width of PML region (in mm)

    // SigmaMaximum, using polynomial grading (Nikolova part 4, p.30)
    electricalConductivityMaximum = -log(reflectionCoefficient0) * (gradingOrder + 1.0) * electricalPermittivity0 * speedOfLight / (2.0 * boundaryWidth);

    // boundaryFactor comes from the polynomial grading equation: sigma_x = sigmaxMaximum * (x/d)^m, where d=width of PML, m=gradingOrder, sigmaxMaximum = electricalConductivityMaximum    (Nikolova part4, p.28)
    //  IMPORTANT: The conductivity (sigma) must use the "average" value at each mesh point as follows:
    //  sigma_x = sigma_Maximum/dx * Integral_from_x0_to_x1 of (x/d)^m dx,  where x0=currentx-0.5, x1=currentx+0.5   (Nikolova part 4, p.32)
    //  integrating gives: sigma_x = (sigmaMaximum / (dx * d^m * m+1)) * ( x1^(m+1) - x0^(m+1) )     (Nikolova part 4, p.32)
    //  the first part is "boundaryFactor", so, sigma_x = boundaryFactor * ( x1^(m+1) - x0^(m+1) )   (Nikolova part 4, p.32)
    // note: it's not exactly clear what the term eps[0] is for. It's probably to cover the case in which eps[0] is not equal to one (ie the main grid area next to the pml boundary is not vacuum)
    boundaryFactor = mediaPermittivity[0] * electricalConductivityMaximum / ( dx * (pow(boundaryWidth,gradingOrder)) * (gradingOrder + 1));

    // build the gradient
    //  caution: if the gradient is built improperly, the PML will not function correctly
    for (i = 0, x = 0.0; i < abcSize; i++, x++) {
        // 0=border between pml and vacuum
        // even: for ex and ey
        x1 = (x + 0.5) * dx;       // upper bounds for point i
        x2 = (x - 0.5) * dx;       // lower bounds for point i
        if (i == 0) {
            gradientConductivity = boundaryFactor * (pow(x1,(gradingOrder+1))  );   //   polynomial grading  (special case: on the edge, 1/2 = pml, 1/2 = vacuum)
        } /* if */
        else {
            gradientConductivity = boundaryFactor * (pow(x1,(gradingOrder+1)) - pow(x2,(gradingOrder+1)) );   //   polynomial grading
        } /* else */
        gradientCa1[i] = exp (-gradientConductivity * dt / (electricalPermittivity0 * mediaPermittivity[0]) );     // exponential time step, Taflove1995 p.77,78
        gradientCb1[i] = (1.0 - gradientCa1[i]) / (gradientConductivity * dx);                                     // ditto, but note sign change from Taflove1995

        // odd: for hzx and hzy
        x1 = (x + 1.0) * dx;       // upper bounds for point i
        x2 = (x + 0.0) * dx;       // lower bounds for point i
        gradientConductivity = boundaryFactor * (pow(x1,(gradingOrder+1)) - pow(x2,(gradingOrder+1)) );   //   polynomial grading
        gradientResistivity = gradientConductivity * (magneticPermeability0 / (electricalPermittivity0 * mediaPermittivity[0]) );  // Taflove1995 p.182  (for no reflection: sigmaM = sigmaE * mu0/eps0)
        gradientDa1[i] = exp(-gradientResistivity * dt / magneticPermeability0);                                   // exponential time step, Taflove1995 p.77,78
        gradientDb1[i] = (1.0 - gradientDa1[i]) / (gradientResistivity * dx);                                      // ditto, but note sign change from Taflove1995
    } /* iForLoop */

    // ex --- front/back
    for (j = 0; j < abcSize; j++) {                            // do coefficients for ex
        for (i = 0; i < xSize; i++) {
            // do coefficients for ex for front and back regions
            caex[i][abcSize - j]         = gradientCa1[j];        // front
            cbex[i][abcSize - j]         = gradientCb1[j];
            caex[i][ySize - abcSize + j] = gradientCa1[j];        // back
            cbex[i][ySize - abcSize + j] = gradientCb1[j];
        } /* iForLoop */
    } /* jForLoop */

    // ey & hzx --- left/right
    for (i = 0; i < abcSize; i++) {                            // do coefficients for ey and hzx
        for (j = 0; j < ySize; j++) {
            // do coefficients for ey for left and right regions
            caey[abcSize - i][j]         = gradientCa1[i];     // left
            cbey[abcSize - i][j]         = gradientCb1[i];
            caey[xSize - abcSize + i][j] = gradientCa1[i];     // right
            cbey[xSize - abcSize + i][j] = gradientCb1[i];
            dahz[abcSize - i - 1][j]     = gradientDa1[i];     // dahz holds dahzx , left, (note that the index is offset by 1 from caey)
            dbhz[abcSize - i - 1][j]     = gradientDb1[i];     // dbhz holds dbhzx         ( ditto )
            dahz[xSize - abcSize + i][j] = gradientDa1[i];     // dahz holds dahzx , right
            dbhz[xSize - abcSize + i][j] = gradientDb1[i];     // dbhz holds dbhzx
        } /* iForLoop */
    } /* jForLoop */

    boundaryIndex = 0;
    // ALERT: special case for hzy:
    // dahzy and dbhzy arrays must be initialized in the same order that they will be accessed in the main time-stepping loop.
    // Therefore it is critical to increment boundaryIndex in the right order for the front/back regions
    // For the left and right regions dahzy,dbhzy use vacuum values, which they are initialized to by default.
    for (i = regionData[1].xStart; i < regionData[1].xStop; i++) {       // front
        for (j = regionData[1].yStart, k=(abcSize - 1); j < regionData[1].yStop; j++, k--) {
            dahzy[boundaryIndex] = gradientDa1[k];
            dbhzy[boundaryIndex] = gradientDb1[k];
            boundaryIndex++;
        } /* jForLoop */
    } /* iForLoop */
    for (i = regionData[2].xStart; i < regionData[2].xStop; i++) {       // back
        for (j = regionData[2].yStart, k=0; j < regionData[2].yStop; j++, k++) {
            dahzy[boundaryIndex] = gradientDa1[k];
            dbhzy[boundaryIndex] = gradientDb1[k];
            boundaryIndex++;
        } /* jForLoop */
    } /* iForLoop */

    // all done with Initialization!

    return;

}
