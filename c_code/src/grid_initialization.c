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
double getMatPlasma(int metalChoice) {
  if ( metalChoice >= 0 ) {
    return materialData[metalChoice][4];
  } else {
    return HUGE_VAL;
  } /* if Block */
}
double getMatDamping(int metalChoice) {
  if ( metalChoice >= 0 ) {
    return materialData[metalChoice][5];
  } else {
    return 0.0;
  } /* if Block */
}
double getMatPermittivity(int metalChoice, double objectIndex) {
  if ( metalChoice >= 0 ) {
    return materialData[metalChoice][0];
  } else {
    return objectIndex*objectIndex;
  } /* if Block */
}
double getMatConductivity(int metalChoice) {
  if ( metalChoice >= 0 ) {
    return materialData[metalChoice][1];
  } else {
    return 0.0;
  } /* if Block */
}
double getMatPermeability(int metalChoice) {
  if ( metalChoice >= 0 ) {
    return materialData[metalChoice][2];
  } else {
    return 1.0;
  }
}
double getMatResistivity(int metalChoice) {
  if ( metalChoice >= 0 ) {
    return materialData[metalChoice][3];
  } else {
    return 0.0;
  } /* if Block */
}

// Function to get the index we're using to reference our "empty" runs by
// this is, unfortunately, the index of the (refractive) index.
int getIndexIndex(double environmentIndex) {
  if (environmentIndex > 4.0) {
    environmentIndex = 4.0;
  } else if (environmentIndex < 1.0) {
    environmentIndex = 1.0;
  } /* if Block */

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


    /* Added for Drude metals so we can treat both the vacuum and object as "Drude"
       materials */
    double  mediaPlasma[MEDIACONSTANT] = {HUGE_VAL, getMatPlasma(metalChoice)}; // Plasma frequency
    double  mediaDamping[MEDIACONSTANT] = {0.0, getMatDamping(metalChoice)}; // Damping constant
    /* End of Drude Metal Addition */

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
    int  i,j,k, boundaryDataSize, media, boundaryIndex,xSizeMain,ySizeMain,numFreqs;
    int  abcSize ;
    double  cylinderDiameter, cylinderRadius, temporaryi,temporaryj,distance2 ;
    int  xCenter,yCenter;
    double  x,x1,x2;
    double  electricalConductivityMaximum, boundaryWidth, gradientConductivity, gradientResistivity, boundaryFactor;
    double  gradientCa1[ABCSIZECONSTANT];
    double  gradientCb1[ABCSIZECONSTANT];
    double  gradientDa1[ABCSIZECONSTANT];
    double  gradientDb1[ABCSIZECONSTANT];
    double  rtau, tau, delay;
    // char  ch;

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
    hzy = AllocateMemory1D(boundaryDataSize, 0.0 );        // for the split-field data for hz in the pml regions

    e2Field = AllocateMemory(xSize, ySize, 0.0);           // E^2 for plotting
    edgeMat = AllocateMemory(xSize, ySize, 0.0);

    /* Polarization Current Fields */
    jx = AllocateMemory(xSize, ySize, 0.0);
    jy = AllocateMemory(xSize, ySize, 0.0);
    exOld = AllocateMemory(xSize, ySize + 1, 0.0);
    eyOld = AllocateMemory(xSize + 1, ySize, 0.0);

    /* Pre-compute values that will be placed in to cjj and cje */
    double cjjTemp[media], cjeTemp[media];
    double nDamping;

    for (i = 0; i < media; i++) { // Schneider 10.57, 10.58, and 10.54 for Ng conversion
        nDamping = 1.0 / (mediaDamping[i] * dt);
        temporary = 1.0 / (2.0 * nDamping);
        cjjTemp[i] = (1.0 - temporary) / (1.0 + temporary);
        // Extra dx on Cje is due to 10.58 having been multiplied through by dx to cancel it out when expressed in it's more traditional form
        cjeTemp[i] = dx * (1.0 / (1.0 + temporary)) * \
         (2.0*courantS*M_PI*M_PI) / (electricalImpedance0 * pow(mediaPlasma[i],2.0));
    } /* iForLoop */

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

    for (i = 0; i < media; i++) {

        /* Original media coefficients
        temporary  = dt * mediaConductivity[i] / (2.0 * electricalPermittivity0 * mediaPermittivity[i] );    // Taflove1995 p.67
        mediaCa[i] = (1.0 - temporary) / (1.0 + temporary);                                                  // ditto
        mediaCb[i] = dt / (electricalPermittivity0 * mediaPermittivity[i] * dx * (1.0 + temporary));         // ditto
        */

        /* Drude material coefficients */
        temp1 = mediaConductivity[i] * dt / (2.0 * electricalPermittivity0 * mediaPermittivity[i]);
        temp2 = cjeTemp[i] * electricalImpedance0 * courantS / (2.0  * mediaPermittivity[i]);
        temp3 = electricalImpedance0 * courantS / mediaPermittivity[i];
        mediaCa[i] = (1.0 - temp1 - temp2) / (1.0 + temp1 + temp2);                                          // See Schneider 10.62
        mediaCb[i] = temp3 / (1.0 + temp1 + temp2);                                                          // ""

        /* Original H updates
        temporary  = dt *  mediaResistivity[i] / (2.0 * magneticPermeability0 * mediaPermeability[i]);       // Taflove1995
        mediaDa[i] = (1.0 - temporary) / (1.0 + temporary);                                                  // ditto
        mediaDb[i] = dt / (magneticPermeability0 * mediaPermeability[i] * dx * (1.0 + temporary));           // ditto
        */
        //temporary = mediaConductivity[i] * dt / (2.0 * electricalPermittivity0 * mediaPermittivity[i]); // Schneider 8.13 - 8.18
        mediaDa[i] = 1.0; // assuming magnetic conductivity is 0
        mediaDb[i] = dt / (dx * magneticPermeability0 * mediaPermeability[i]);
    } /* iForLoop */


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

    caex = AllocateMemory(xSize, ySize, mediaCa[0] );     // note: don't need to allocate for pec region, as it is not evaluated
    cbex = AllocateMemory(xSize, ySize, mediaCb[0] );     // also: Initialize the entire grid to vacuum.
    caey = AllocateMemory(xSize, ySize, mediaCa[0] );
    cbey = AllocateMemory(xSize, ySize, mediaCb[0] );
    dahz = AllocateMemory(xSize, ySize, mediaDa[0] );
    dbhz = AllocateMemory(xSize, ySize, mediaDb[0] );
    dahzy = AllocateMemory1D(boundaryDataSize, mediaDa[0] );        // for the split-field data for hz in the pml regions
    dbhzy = AllocateMemory1D(boundaryDataSize, mediaDb[0] );        // for the split-field data for hz in the pml regions

    /* Initialize polarization current */
    cjj = AllocateMemory(xSize, ySize, 1.0);  // Initialize to 1 for free space
    cje = AllocateMemory(xSize, ySize, 0.0);  // Initialize to 0 for free space

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

    // Add structure:
    for(i = 0; i < xSize; i++) {
      for (j = 0; j < ySize; j++) {
        if (object_locs[i][j] > 0.5) {
          // Material Constants:
          caex[i][j] = mediaCa[1];
          cbex[i][j] = mediaCb[1];
          caey[i][j] = mediaCa[1];
          cbey[i][j] = mediaCb[1];

          /* Polarization Current Constants: */
          cjj[i][j] = cjjTemp[1];
          cje[i][j] = cjeTemp[1];
        }
      } /* jForLoop */
    } /* iForLoop */
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
    regionData[1].xStart = 0;                          // front grid
    regionData[1].xStop  = xSize;
    regionData[1].yStart = 0;
    regionData[1].yStop  = abcSize;
    regionData[2].xStart = 0;                          // back grid
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
