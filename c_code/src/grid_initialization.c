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
    double  drudePlasma[MEDIACONSTANT] = {0.0};
    double  drudeDamping[MEDIACONSTANT] = {HUGE_VAL}; // initialize vacuum to have huge damping
    complex double  mediaA[MEDIACONSTANT][MAX_POLES] = {0.0};
    complex double  mediaC[MEDIACONSTANT][MAX_POLES] = {0.0};
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

    courantS = 1.0/2.0;
    dt = courantS * dx / speedOfLight;

    absConst = (speedOfLight*dt - dx) / (speedOfLight*dt + dx);
//printf( "dx: %f\n", dx );
//printf( "dt: %f\n", dt );
//printf( "courantS: %f\n", courantS );
    /***********************************************************************/
    //     Grid parameters
    /***********************************************************************/

    /*xSizeMain = 300;                              // number of main grid cells in x-direction
    ySizeMain = 250;*/
    xSizeMain = 200;
    ySizeMain = 150;                               // number of main grid cells in y-direction
    abcSize = ABCSIZECONSTANT;                    // thickness of PML region
    xSize = xSizeMain + 2 * abcSize;              // number of total grid cells in x-direction
    ySize = ySizeMain + 2 * abcSize;              // number of total grid cells in y-direction

    boundaryDataSize  = 2 * xSize * abcSize;                      // front edge + back edge
    boundaryDataSize += 2 * (abcSize * (ySize - 2 * abcSize));    // left + right edges

    xSource = 50 + abcSize;                          //location of z-directed hard source
    //ySource = 50 + abcSize;                          //location of z-directed hard source
    //xSource = 20 + abcSize;                       //location of z-directed hard source
    ySource = ySize / 2;                          //location of z-directed hard source

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

    heConst = AllocateMemory(xSize, ySize, dt/(dx*magneticPermeability0));
    ehConst = AllocateMemory(xSize, ySize, dt/(dx*electricalPermittivity0));
    eqConst = AllocateMemory(xSize, ySize, 4.0*dt/electricalPermittivity0);
    qxSum   = AllocateMemory(xSize, ySize, 0.0);
    qySum   = AllocateMemory(xSize, ySize, 0.0);
    qConst1 = AllocateComplexMemory3D(number_poles, xSize, ySize, initArray);
    qConst2 = AllocateComplexMemory3D(number_poles, xSize, ySize, initArray);
    qSumC   = AllocateComplexMemory3D(number_poles, xSize, ySize, initArray);
    ABConst = AllocateMemory(xSize, ySize, dt*dt/(4.0*magneticPermeability0*electricalPermittivity0*dx*dx));


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
          ABConst[i][j] /= mediaPermittivity[1]*mediaPermeability[1];
          ehConst[i][j] /= mediaPermittivity[1];
          eqConst[i][j] /= mediaPermittivity[1];
          heConst[i][j] /= mediaPermeability[1];

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


    /***********************************************************************/
    //     Fill the PML with lossy terms
    /***********************************************************************/

    // The most important part of the PML fdtd simulation is getting the
    // PML Coefficients correct. Which requires getting the correct PML gradient and
    // positioning the coefficients correctly on the x-y grid.

    // ALERT: It is possible to make a mistake here, and yet the simulation may appear to
    // be working properly. However a detailed analysis of reflections off the PML
    // will show they may be (much) larger than those for a correctly designed PML.

    gradingOrder = 4.0;                             // for PML, (m) was 2;  optimal values: 2 <= m <= 6,  Nikolova part4 p.29

    double gradKx[ABCSIZECONSTANT] = {0.0};
    double gradKy[ABCSIZECONSTANT] = {0.0};
    double kMax = 4.0;
    PMLkx = AllocateMemory(xSize, ySize, 1.0); // Multipled anywhere we have d^2 / dx^2 derivatives
    PMLky = AllocateMemory(xSize, ySize, 1.0); // " but for y-derivatives
    //PMLkz = AllocateMemory(xSize, ySize, 1.0); // This is always 1.0, but it's nice to have for symmetry


    boundaryWidth = (double  )abcSize * dx;    // width of PML region (in mm)
    double prefactor = (kMax-1.0)/(dx * pow(boundaryWidth,gradingOrder)) * (gradingOrder + 1.0);

    // build the gradient
    //  caution: if the gradient is built improperly, the PML will not function correctly
    for (i = 0, x = 0.0; i < abcSize; i++, x++) {
        x1 = (x + 0.5) * dx;       // upper bounds for point i
        x2 = (x - 0.5) * dx;       // lower bounds for point i
        if (i == 0) {
            temporary = 1.0 + prefactor * (pow(x1,(gradingOrder+1))  );   //   polynomial grading  (special case: on the edge, 1/2 = pml, 1/2 = vacuum)
        } /* if */
        else {
            temporary = 1.0 + prefactor * (pow(x1,(gradingOrder+1)) - pow(x2,(gradingOrder+1)) );   //   polynomial grading
        } /* else */
        gradKy[i] = temporary;

        x1 = (x + 1.0) * dx;       // upper bounds for point i
        x2 = (x + 0.0) * dx;       // lower bounds for point i
        gradKx[i] = 1.0 + prefactor * (pow(x1,(gradingOrder+1)) - pow(x2,(gradingOrder+1)) );   //   polynomial grading

    } /* iForLoop */

    printf("Order: %f\n", gradingOrder);
    printf("gradKx[0]: %.17g\n", gradKx[0]);
    printf("gradKx[7]: %.17g\n", gradKx[7]);

    // Now, loop through PML regions and set PML-k terms where appropriate:
    // Left/Right PML Region:
    /*for (i = 0; i < abcSize; i++) {
      for (j = 0; j < ySize; j++) {
        PMLkx[abcSize - 1 - i][j]     = 1.0/gradKx[i];
        PMLkx[xSize - abcSize + i][j] = 1.0/gradKx[i];
      }
    }

    // Top/Bottom
    for (j = 0; j < abcSize; j++) {
      for (i = 0; i < xSize; i++) {
        PMLky[i][abcSize - 1 - j]     = 1.0/gradKy[j];
        PMLky[i][ySize - abcSize + j] = 1.0/gradKy[j];
      }
    }*/
    /*
    printf("PMLkx[7][133]: %.17g\n", PMLkx[7][133]);
    printf("ABCConst*PMLkx[7][133]: %.17g\n", ABConst[7][133]*PMLkx[7][133]);
    printf("PMLkx[0][0]: %.17g\n", PMLkx[0][0]);
    printf("ABCConst*PMLkx[0][0]: %.17g\n", ABConst[0][0]*PMLkx[0][0]);
    */

    // Finally, set values for our tridiagonal solver(s)
    // This needs to be done last as they take the values from the PML where necessary

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
        ahz[i][j] = -1.0 * PMLkx[i][j]*ABConst[i][j];
        bhz[i][j] = iConst1[i][j] + 2.0*PMLkx[i][j]*ABConst[i][j];
        chz[i][j] = ahz[i][j];

        aex[i][j] = -1.0 * PMLky[i][j]*ABConst[i][j];
        bex[i][j] = iConst1[i][j] + 2.0*PMLky[i][j]*ABConst[i][j];
        cex[i][j] = aex[i][j];
      }
    }

    // Add Mur ABC conditions from "Mur Absorbing Boundary Condition for 2-D
    // Leapgfrog ADI-FDTD Method", Gan and Tan, IEEE Asia-Pacific Conf. 2012

    // As Ex depends on j-terms:
    for (j = 0; j < ySize; j++) {
      aex[0][j] = 0.0;
      bex[0][j] = 1.0;
      cex[0][j] = -1.0*absConst;

      aex[xSize-1][j] = -1.0*absConst;
      bex[xSize-1][j] = 1.0;
      cex[xSize-1][j] = 0.0;
    }

    // While Hz depends on i-terms:
    for (i = 0; i < xSize; i++) {
      ahz[i][0] = 0.0;
      bhz[i][0] = 1.0;
      chz[i][0] = -1.0*absConst;

      ahz[i][ySize-1] = -1.0*absConst;
      bhz[i][ySize-1] = 1.0;
      chz[i][ySize-1] = 0.0;
    }

    printf("ahz: %.17g\n", ahz[4][0]);
    printf("bhz: %.17g\n", bhz[4][0]);

    printf("aex: %.17g\n", aex[4][0]);
    printf("bex: %.17g\n", bex[4][0]);


    // all done with Initialization!

    return;

}
