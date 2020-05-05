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
  double objectSize, double environmentIndex, double objectIndex )
{

    double  mediaPermittivity[MEDIACONSTANT] = {environmentIndex*environmentIndex, getMatPermittivity(metalChoice, objectIndex)};    // eps, index=0 is for vacuum, index=1 is for the metallic cylinder
    double  mediaConductivity[MEDIACONSTANT] = {0.0, getMatConductivity(metalChoice)}; // sig,
    double  mediaPermeability[MEDIACONSTANT] = {1.0, getMatPermeability(metalChoice)};    // mur
    double  mediaResistivity[MEDIACONSTANT] = {0.0, getMatResistivity(metalChoice)};     // sim
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
    double  reStretched[ABCSIZECONSTANT] = {0.0}; // Re part of stretched coordinate system in PML (eq. 34 Prokopidis and Zografopoulos)
    double  imStretched[ABCSIZECONSTANT] = {0.0};
    double  kappaPML[ABCSIZECONSTANT] = {0.0};
    double  sigmaPML[ABCSIZECONSTANT] = {0.0};

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
    xSource = 15 + abcSize;                       //location of z-directed hard source
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

    // Number of poles in our dielectric function
    number_poles = materialData[metalChoice].num_poles;

    // Constants for our approximation:
    double **Ap,**Omegap,**phip,**Gammap;
    Ap = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);
    Omegap = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);
    phip = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);
    Gammap = AllocateMemory(MEDIACONSTANT, number_poles, 0.0);

    for (p = 0; p < number_poles; p++) {
      Ap[1][p] = materialData[metalChoice].params[p].bigA;
      Omegap[1][p] = materialData[metalChoice].params[p].Omega;
      phip[1][p] = materialData[metalChoice].params[p].phi;
      Gammap[1][p] = materialData[metalChoice].params[p].Gamma;
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
    

    e2Field = AllocateMemory(xSize, ySize, 0.0);           // E^2 for plotting
    edgeMat = AllocateMemory(xSize, ySize, 0.0);

    /* Polarization Current Fields */
    double initArray[12] = {0.0};
    px = AllocateMemory3D(number_poles, xSize, ySize, initArray);
    py = AllocateMemory3D(number_poles, xSize, ySize, initArray);
    pxOld = AllocateMemory3D(number_poles, xSize, ySize, initArray);
    pyOld = AllocateMemory3D(number_poles, xSize, ySize, initArray);
    exOld = AllocateMemory(xSize, ySize + 1, 0.0);
    exOld2 = AllocateMemory(xSize, ySize + 1, 0.0);
    eyOld = AllocateMemory(xSize + 1, ySize, 0.0);
    eyOld2 = AllocateMemory(xSize + 1, ySize, 0.0);

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
    double AZero,AOne,BZero,BOne,BTwo,cP;
    double nDamping;

    // Magnetic field update constants:
    for (i = 0; i < media; i++) {
        //temporary = mediaConductivity[i] * dt / (2.0 * electricalPermittivity0 * mediaPermittivity[i]); // Schneider 8.13 - 8.18
        mediaDa[i] = 1.0; // assuming magnetic conductivity is 0
        mediaDb[i] = dt / (dx * magneticPermeability0 * mediaPermeability[i]);
    } /* iForLoop */

    /* For terms and formulation, see: A Unified FDTD/PML Scheme Based on Critical
      Points for Accurate Studies of Plasmonic Structures. Prokopidis and
      Zografopoulos, JOURNAL OF LIGHTWAVE TECHNOLOGY, VOL. 31, NO. 15, AUGUST 1, 2013 */

    // Electric / Polarization field update constants:
    // Critical Points terms:
    for (i = 0; i < media; i++) {
      for (p = 0; p < number_poles; p++){
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
        c4TempSum[i] += c4[i][p];
        c5TempSum[i] += c5[i][p];
      } /* pForLoop */
      // See eq. 23 in the paper:
      c3TempSum[i] += electricalPermittivity0*mediaPermittivity[i];
      c4TempSum[i] -= electricalPermittivity0*mediaPermittivity[i];
    } /* iForLoop */

    // Drude Terms (assuming one Drude term)
    /*for (i = 0; i < media; i++) {
      // Drude Points Terms
      dQ = 1.0/(dt*dt) + gamma/(2.0*dt);
      DOne = 2.0/(dQ*dt*dt);
      DTwo = (gamma/(2.0*dt) - 1.0/(dt*dt))/dQ;
      DThree = electricalPermittivity0*mediaPlasma[i]*mediaPlasma[i]/(4.0*dQ);
      DFour = DThree*2.0;
      DFive = DThree;
    }*/ /* iForLoop */

    /***********************************************************************/
    //     Grid Coefficients
    /***********************************************************************/
    /*printf("mediaCa[0]: %f\n", mediaCa[0]);
    printf("mediaCa[1]: %f\n", mediaCa[1]);
    printf("mediaCb[0]: %f\n", mediaCb[0]);
    printf("mediaCb[1]: %f\n", mediaCb[1]);
    printf("mediaDa[0]: %f\n", mediaDa[0]);
    printf("mediaDa[1]: %f\n", mediaDa[1]);
    printf("mediaDb[0]: %f\n", mediaDb[0]);
    printf("mediaDb[1]: %f\n", mediaDb[1]);*/
    //     Initialize entire grid to free space

    // Polarization grid values:
    c1SumX = AllocateMemory(xSize, ySize, c3TempSum[0]); // No temp sums for c1 and c2, and this is 0 anyway
    c2SumX = AllocateMemory(xSize, ySize, c3TempSum[0]);
    c1SumY = AllocateMemory(xSize, ySize, c3TempSum[0]);
    c2SumY = AllocateMemory(xSize, ySize, c3TempSum[0]);
    c3Sum = AllocateMemory(xSize, ySize, c3TempSum[0]);
    c4Sum = AllocateMemory(xSize, ySize, c4TempSum[0]);
    c5Sum = AllocateMemory(xSize, ySize, c5TempSum[0]);


    // Make so we don't have to loop over each array twice on initialization...
    c1Grid = AllocateMemory3D(number_poles, xSize, ySize, c1[0]);
    c2Grid = AllocateMemory3D(number_poles, xSize, ySize, c2[0]);
    c3Grid = AllocateMemory3D(number_poles, xSize, ySize, c3[0]);
    c4Grid = AllocateMemory3D(number_poles, xSize, ySize, c4[0]);
    c5Grid = AllocateMemory3D(number_poles, xSize, ySize, c5[0]);

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


    eGrad1 = AllocateMemory1D(boundaryDataSize, 0.0 );        // PML data
    hGrad1 = AllocateMemory1D(boundaryDataSize, 0.0 );
    eGrad2 = AllocateMemory1D(boundaryDataSize, 0.0 );
    hGrad2 = AllocateMemory1D(boundaryDataSize, 0.0 );
    eGrad3 = AllocateMemory1D(boundaryDataSize, 0.0 );
    hGrad3 = AllocateMemory1D(boundaryDataSize, 0.0 );

    pmlSx = AllocateMemory1D(boundaryDataSize, 0.0 );
    pmlSy = AllocateMemory1D(boundaryDataSize, 0.0 );
    pmlTz = AllocateMemory1D(boundaryDataSize, 0.0 );
    pmlSxOld = AllocateMemory1D(boundaryDataSize, 0.0 );
    pmlSyOld = AllocateMemory1D(boundaryDataSize, 0.0 );
    pmlTzOld = AllocateMemory1D(boundaryDataSize, 0.0 );

    rx = AllocateMemory1D(boundaryDataSize, 0.0 );
    rxOld = AllocateMemory1D(boundaryDataSize, 0.0 );
    rxOld2 = AllocateMemory1D(boundaryDataSize, 0.0 );
    ry = AllocateMemory1D(boundaryDataSize, 0.0 );
    ryOld = AllocateMemory1D(boundaryDataSize, 0.0 );
    ryOld2 = AllocateMemory1D(boundaryDataSize, 0.0 );
    bz = AllocateMemory1D(boundaryDataSize, 0.0 );
    /* Values for H-field updates: */
    dahz = AllocateMemory(xSize, ySize, mediaDa[0] );
    dbhz = AllocateMemory(xSize, ySize, mediaDb[0] );

    /* Array to track where our object is/is not */
    object_locs = AllocateMemory(xSize, ySize, 0.0); // This really shouldn't be an array of doubles -- we're wasting memory here

    // Initialize structure functions:
    xCenter = xSize - abcSize - 30; // In grid units
    yCenter = ySize / 2; // ""
    structInit(xCenter, yCenter);
printf("Strucutre Init...\n" );
    // Switch Block to pick structure geometry (default is no object):
    switch (objectChoice) {
      case 0: // Disk
        addDisk(g, objectSize * dx/(2.0 * dxnm)); // 10 * dx radius disk
        break;

      case 1: // Block
        addRect(g, 20.0 * dx, objectSize * dx/dxnm);
        break;

      case 2: // Triangle
        addTriangle(g, objectSize * dx/dxnm);
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

          dahz[i][j] = mediaDa[1];
          dbhz[i][j] = mediaDb[1];

          for(p = 0; p < number_poles; p++) {
            /* Polarization Constants: */
            c1Grid[p][i][j] = c1[1][p];
            c2Grid[p][i][j] = c2[1][p];
            c3Grid[p][i][j] = c3[1][p];
            c4Grid[p][i][j] = c4[1][p];
            c5Grid[p][i][j] = c5[1][p];
          } /* pForLoop */
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

    boundaryWidth = (double  )abcSize;    // width of PML region

    // SigmaMaximum, using polynomial grading (Nikolova part 4, p.30)
    electricalConductivityMaximum = -log(reflectionCoefficient0) * (gradingOrder + 1.0) * electricalPermittivity0 * speedOfLight / (2.0 * boundaryWidth);

    kappaMaximum = 2.0; // ????? I don't know how to pick this.

    alphaPML = 2.0; // ????? Also don't know how to pick this. Likely tuning...

    for (i = 0; i < ABCSIZECONSTANT; i++) {
      temporary = (double  )i;
      sigmaPML[i] = electricalConductivityMaximum * pow((temporary/boundaryWidth), gradingOrder);
      kappaPML[i] = 1.0 + (kappaMaximum - 1.0) * pow((temporary/boundaryWidth), gradingOrder);
    } /* i forLoop */

    boundaryIndex = 0;
    // Bottom:
    for (i = regionData[1].xStart; i < regionData[1].xStop; i++) {
      for (j = regionData[1].yStart, k=(abcSize - 1); j < regionData[1].yStop; j++, k--) {
        denominator = electricalPermittivity0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt;
        eGrad1[boundaryIndex] = ( electricalPermittivity0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt ) / denominator;
        eGrad2[boundaryIndex] = ( electricalPermittivity0 + 0.5*alphaPML*dt )/denominator;
        eGrad3[boundaryIndex] = ( electricalPermittivity0 - 0.5*alphaPML*dt )/denominator;

        denominator = magneticPermeability0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt;
        hGrad1[boundaryIndex] = ( magneticPermeability0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt ) / denominator;
        hGrad2[boundaryIndex] = ( magneticPermeability0 + 0.5*alphaPML*dt )/denominator;
        hGrad3[boundaryIndex] = ( magneticPermeability0 - 0.5*alphaPML*dt )/denominator;
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */

    // Top:
    for (i = regionData[2].xStart; i < regionData[2].xStop; i++) {
      for (j = regionData[2].yStart, k=0; j < regionData[2].yStop; j++, k++) {
        denominator = electricalPermittivity0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt;
        eGrad1[boundaryIndex] = ( electricalPermittivity0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt ) / denominator;
        eGrad2[boundaryIndex] = ( electricalPermittivity0 + 0.5*alphaPML*dt )/denominator;
        eGrad3[boundaryIndex] = ( electricalPermittivity0 - 0.5*alphaPML*dt )/denominator;

        denominator = magneticPermeability0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt;
        hGrad1[boundaryIndex] = ( magneticPermeability0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt ) / denominator;
        hGrad2[boundaryIndex] = ( magneticPermeability0 + 0.5*alphaPML*dt )/denominator;
        hGrad3[boundaryIndex] = ( magneticPermeability0 - 0.5*alphaPML*dt )/denominator;
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */

    // Left:
    for (i = regionData[3].xStart, k=(abcSize - 1); i < regionData[3].xStop; i++, k--) {
      for (j = regionData[3].yStart; j < regionData[3].yStop; j++) {
        denominator = electricalPermittivity0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt;
        eGrad1[boundaryIndex] = ( electricalPermittivity0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt ) / denominator;
        eGrad2[boundaryIndex] = ( electricalPermittivity0 + 0.5*alphaPML*dt )/denominator;
        eGrad3[boundaryIndex] = ( electricalPermittivity0 - 0.5*alphaPML*dt )/denominator;

        denominator = magneticPermeability0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt;
        hGrad1[boundaryIndex] = ( magneticPermeability0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt ) / denominator;
        hGrad2[boundaryIndex] = ( magneticPermeability0 + 0.5*alphaPML*dt )/denominator;
        hGrad3[boundaryIndex] = ( magneticPermeability0 - 0.5*alphaPML*dt )/denominator;
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */

    // Right:
    for (i = regionData[4].xStart, k=0; i < regionData[4].xStop; i++, k++) {
      for (j = regionData[4].yStart; j < regionData[4].yStop; j++) {
        denominator = electricalPermittivity0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt;
        eGrad1[boundaryIndex] = ( electricalPermittivity0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt ) / denominator;
        eGrad2[boundaryIndex] = ( electricalPermittivity0 + 0.5*alphaPML*dt )/denominator;
        eGrad3[boundaryIndex] = ( electricalPermittivity0 - 0.5*alphaPML*dt )/denominator;

        denominator = magneticPermeability0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt;
        hGrad1[boundaryIndex] = ( magneticPermeability0*kappaPML[k] + \
          0.5*(sigmaPML[i] + alphaPML*kappaPML[i])*dt ) / denominator;
        hGrad2[boundaryIndex] = ( magneticPermeability0 + 0.5*alphaPML*dt )/denominator;
        hGrad3[boundaryIndex] = ( magneticPermeability0 - 0.5*alphaPML*dt )/denominator;
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */

    // all done with Initialization!

    return;

}
