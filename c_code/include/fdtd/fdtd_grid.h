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

#ifndef FDTD_GRID
#define FDTD_GRID

#include <complex.h>

#define  ABCSIZECONSTANT    (25)                      // thickness of PML region
#define  MEDIACONSTANT    (2)                        // number of different media, ie 2: vacuum, metal object
#define  NUMBEROFITERATIONCONSTANT    (1500)          // Number of timesteps
#define  NUMBEROFREGIONS    (5)                      // center(main), front, back, left, right

#define  NUMBERDFTFREQS    (50)                      // Number of frequencies to compute DFTs at
#define  DFTPADDEDTIME    (16384)

struct RegionDataValues{
  int xStart,yStart,xStop,yStop;
};

//typedef struct RegionDataValues RegionDataValues;

struct Grid {
  // These global variables are used by ExecuteFdtd()

  int  xSize ;                           // number of total grid cells in x-direction  (left+main+right)
  int  ySize ;                           // number of total grid cells in y-direction  (front+main+back)
  int  maximumIteration;                 // how long to run the simulation
  int  xSource, ySource;                 // location of z-directed hard source

  double dx,dt,courantS,speedOfLight,pi;
  double envIndex; // Refractive index of the environment

  double  **ex;      // the fields
  double  **ey;      //  ""
  double  **hz;      //  ""

  complex double  ***qx;
  complex double  ***qy;
  complex double  ***qSumC;
  int number_poles;
  double **heConst;
  double **ehConst;
  double **eqConst;
  double **ABConst;
  complex double ***qConst1;
  complex double ***qConst2;
  double **qxSum;
  double **qySum;
  double **iConst1;
  double **iConst2;

  // Tridiagonal values:
  double  **aex;
  double  **bex;
  double  **cex;
  double  **ahz;
  double  **bhz;
  double  **chz;

  double  **exOld;    // Matrix to store old Ex values for Drude metals
  double  **eyOld;    // "" Ey values
  double  **hzOld;

  // ABC update value:
  double  absConst;

  // Values for tracking DFT
  int reflXPos,tranXPos;          // Positions of the line monitors used to find reflected / transmitted fields
  double *reflDFT;                // Values for reflection DFT at some number of wavelengths
  double *tranDFT;                // "" for transmission
  int kList[NUMBERDFTFREQS];       // List storing the wave vectors we are taking the DFTs at
  double wavelengthList[NUMBERDFTFREQS]; // "" wavelengths

  // Arrays to store DFT running sums as a function of position and frequency
  // Need to store Ey and Hz at each point on the line we are integrating over

  double **reEyReflDFT;
  double **imEyReflDFT;
  double **reEyTranDFT;
  double **imEyTranDFT;

  double **reHzReflDFT;
  double **imHzReflDFT;
  double **reHzTranDFT;
  double **imHzTranDFT;

  // E^2 field for plotting:
  double **e2Field;

  // Matrix containing edge locations:
  double **edgeMat;

  // Array to track where our object is/is not:
  double **object_locs;
  int  objectXMax; // Variables to store greatest extent of the object so we don't have to update Q over the entire simulation space
  int  objectYMax;
  int  objectXMin;
  int  objectYMin;

  // Index for our "empty" runs to track what index we're using now:
  int refractiveIndexIndex; // unfortunately named...

  double  sourceValue[NUMBEROFITERATIONCONSTANT];  // holds the pre-calculated values of the source, for the run of the simulation

  struct RegionDataValues  regionData[NUMBEROFREGIONS];   // for calculating Hz (includes information about PML regions)

  // ints to track iterations when using emscripten:
  int timeStep;
  int interval;
};
//typedef struct Grid Grid;



#endif
