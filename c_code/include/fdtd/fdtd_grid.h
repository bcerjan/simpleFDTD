#ifndef FDTD_GRID
#define FDTD_GRID

#define  ABCSIZECONSTANT    (8)                      // thickness of PML region
#define  MEDIACONSTANT    (2)                        // number of different media, ie 2: vacuum, metal object
#define  NUMBEROFITERATIONCONSTANT    (8000)          // Number of timesteps
#define  NUMBEROFREGIONS    (5)                      // center(main), front, back, left, right

#define  NUMBERDFTFREQS    (20)                      // Number of frequencies to compute DFTs at

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

  double  **ex;      // the fields
  double  **ey;      //  ""
  double  **hz;      //  ""
  double  **caex;    // fdtd coefficents
  double  **cbex;    //  ""
  double  **caey;    //  ""
  double  **cbey;    //  ""
  double  **dahz;    // in pml regions this holds dahzx
  double  **dbhz;    //          " "              dbhzx

  double  *hzy;      // for pml split-field abc, note: hzx is derived "on the fly" from hz - hzy
  double  *dahzy;    // pml coefficient
  double  *dbhzy;    //    ""

  // Values for Drude Metals
  double  **jx;       // x-direction polarization current (actually delta*Jx)
  double  **jy;       // "" for Jy
  double  **cjj;      // Matrix containing damping values (10.57 in Schneider)
  double  **cje;      // Matrix containing plasma frequency values (10.58 in Schneider)
  double  **exOld;    // Matrix to store old Ex values for Drude metals
  double  **eyOld;    // "" Ey values

  // Values for tracking DFT
  int reflXPos,tranXPos;          // Positions of the line monitors used to find reflected / transmitted fields
  double *reflDFT;                // Values for reflection DFT at some number of wavelengths
  double *tranDFT;                // "" for transmission
  int kList[NUMBERDFTFREQS];       // List storing the wave vectors we are taking the DFTs at
  double wavelengthList[NUMBERDFTFREQS]; // "" wavelengths

  // Testing storing Re/Im of Ey and Hz instead of alternate running sum...
  double **reEyDFT;
  double **imEyDFT;
  double **reHzDFT;
  double **imHzDFT;

  // E^2 field for plotting:
  double **e2Field;

  // Array to track where our object is/is not:
  double **object_locs;

  double  sourceValue[NUMBEROFITERATIONCONSTANT];  // holds the pre-calculated values of the source, for the run of the simulation

  struct RegionDataValues  regionData[NUMBEROFREGIONS];   // for calculating Hz (includes information about PML regions)

  // ints to track iterations when using emscripten:
  int timeStep;
  int interval;
};
//typedef struct Grid Grid;



#endif
