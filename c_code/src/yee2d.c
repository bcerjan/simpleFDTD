
/* ********************************************************************
    yee2d.c : 2-D FDTD TE code with PML absorbing boundary conditions
   ********************************************************************

   A 2-D FDTD TE implementation of Berengers split-field PML ABC
    version 1.0,  4/4/2008,  Doug Neubauer

   Adapted from the 2-d FDTD TE code of Dr. Susan C. Hagness.

     This program implements the finite-difference time-domain
     solution of Maxwell's curl equations over a two-dimensional
     Cartesian space lattice comprised of uniform square grid cells.

     To illustrate the algorithm, a 6-cm-diameter metal cylindrical
     scatterer in free space is modeled. The source excitation is
     a Gaussian pulse with a carrier frequency of 5 GHz.

     The grid resolution (dx = 3 mm) was chosen to provide 20 samples
     per wavelength at the center frequency of the pulse (which in turn
     provides approximately 10 samples per wavelength at the high end
     of the excitation spectrum, around 10 GHz).

     The computational domain is truncated using the perfectly matched
     layer (PML) absorbing boundary conditions.  The formulation used
     in this code is based on the original split-field Berenger PML. The
     PML regions are labeled as shown in the following diagram:

                                                             xSize-1,ySize-1
                                                            /
            ----------------------------------------------.
           |  |                BACK PML                |  |
            ----------------------------------------------
           |L |                                        | R|
           |E |                                        | I|
           |F |                                        | G|
           |T |                                        | H|
           |  |                MAIN GRID               | T|
           |P |                                        |  |
           |M |                                        | P|
           |L |                                        | M|
           |  |                                        | L|
            ----------------------------------------------
           |. |                FRONT PML               |  |
           /----------------------------------------------
         0,0


   Below: Detailed view of the Yee grid...  (where N=xSize, M=ySize, see below)
   Note: an extra column of ey on right edge and an extra row of ex on back edge
   Note: ey at x=0 and x=N are PEC and ex at y=0 and y=M are PEC.

     (0,M)                                                             (N-1,M)
    ___ex___  ___ex___  ___ex___  ___ex___  ..  ___ex___  ___ex___  ___ex___
    ......... ......... ......... .........     ......... ......... .........
   |        .|        .|        .|        .    |        .|        .|        .|
   | 0,M-1  .|        .|        .|        .    |        .|        .| N-1,M-1.|(N,M-1)
   ey  hz   .ey  hz   .ey  hz   .ey  hz   .    ey  hz   .ey  hz   .ey  hz   .ey
   |        .|        .|        .|        .    |        .|        .|        .|
   |___ex___.|___ex___.|___ex___.|___ex___.    |___ex___.|___ex___.|___ex___.|
    ......... ......... ......... ......... ..  ......... ......... .........
   |        .|        .|        .|        .    |        .|        .|        .|
   |  0,M-2 .|        .|        .|        .    |        .|        .| N-1,M-2.|
   ey  hz   .ey  hz   .ey  hz   .ey  hz   .    ey  hz   .ey  hz   .ey  hz   .ey
   |        .|        .|        .|        .    |        .|        .|        .|
   |___ex___.|___ex___.|___ex___.|___ex___. .. |___ex___.|___ex___.|___ex___.|
   .                                      .    .                             .
   .                                      .    .                             .
    ......... ......... ......... ......... ..  ......... ......... .........
   |        .|        .|        .|        .    |        .|        .|        .|
   |  0,2   .|        .|        .|        .    |        .|        .| N-1,2  .|
   ey  hz   .ey  hz   .ey  hz   .ey  hz   .    ey  hz   .ey  hz   .ey  hz   .ey
   |        .|        .|        .|        .    |        .|        .|        .|
   |___ex___.|___ex___.|___ex___.|___ex___.    |___ex___.|___ex___.|___ex___.|
    ......... ......... ......... ......... ..  ......... ......... .........
   |        .|        .|        .|        .    |        .|        .|        .|
   |  0,1   .|        .|        .|        .    |        .|        .| N-1,1  .|
   ey  hz   .ey  hz   .ey  hz   .ey  hz   .    ey  hz   .ey  hz   .ey  hz   .ey
   |        .|        .|        .|        .    |        .|        .|        .|
   |___ex___.|___ex___.|___ex___.|___ex___.    |___ex___.|___ex___.|___ex___.|
    ......... ......... ......... ......... ..  ......... ......... .........
   |        .|        .|        .|        .    |        .|        .|        .|
   |  0,0   .|  1,0   .|   2,0  .|  3.0   .    |  N-3,0 .|  N-2,0 .| N-1,0  .|(N,0)
   ey  hz   .ey  hz   .ey  hz   .ey  hz   .    ey  hz   .ey  hz   .ey  hz   .ey
   |        .|        .|        .|        .    |        .|        .|        .|
   |___ex___.|___ex___.|___ex___.|___ex___. .. |___ex___.|___ex___.|___ex___.|


     PML Reflection Simulation Results:
     ------------------------------------------
     #layer      r0      gradingOrder     dB
     ------    ------    ------------    -----
       8       1.0e-7        2           -74.4
       8       1.0e-7        3           -81.2      <<== use this one
       8       1.0e-7        4           -81.15
      10       1.0e-7        3           -90.3
      10       1.0e-8        3           -88.7
       6       1.0e-7        3           -71.8
       4       1.0e-7        3           -51.8
     Results are consistant with the version from Dr. Susan C. Hagness.

 In the spirit of the ToyFdtd programs, a point was made to try to heavily comment the source code.

********************************************************************** */

/*
  This version of the code has been modified by Ben Cerjan to make it compatible
  with emscripten to generate a wasm file for embedding in a web page.

  Light is incident in the x-direction from left to right

  InitializeFDTD( int MetalChoice );
    Allowed metals:
      - 0 = Aluminum
      - 1 = Gold
      - 2 = Silver
      - 3 = Copper

  TODO:
    -Allow specification of Metal (and modelling them as Drude Metal, see chap 10 in Schneider)
	= need wp and gamma (damping) for the metals
        = Done! But need to look up real values for relevant materials and store them
    -Allow specification of Geometry
    -Adjust to visible wavelengths
        = Done! (needs testing)
    -Adjust to Ricker source (from Gaussian) and make line instead of point source
    -Add running dft to compute Refl/Trans (and then find Abs)
        = Computing Refl/Trans now
    -Put empty results into a .h file so they can be accessed later
        = Need to do this
    -Add drawing of fields to screen (need to modify drawing function to work
     with memory pointers instead of pass-by-reference)
       = Done!
    -Figure out how to output and plot R/T/A values

  OPTIMIZATIONS:
    -Only track Polarization currents where it is non-zero (i.e. in a metal) instead of everywhere
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fdtd_macro.h"
#include "sdl_funcs.h"

// standard C memory allocation for 2-D array
double  **AllocateMemory (int  imax, int  jmax, double  initialValue)
{
    int  i,j;
    double  **pointer;
    pointer = (double **)malloc(imax * sizeof(double *));
    for (i = 0; i < imax; i++) {
        pointer[i] = (double *)malloc(jmax * sizeof(double));
        for (j = 0; j < jmax; j++) {
            pointer[i][j] = initialValue;
        } /* jForLoop */
    } /* iForLoop */
    return(pointer) ;
}

// standard C memory allocation for 1-D array
double  *AllocateMemory1D (int  size, double  initialValue)
{
    int  j;
    double  *pointer;
    pointer = (double *)malloc(size * sizeof(double));
    for (j = 0; j < size; j++) {
        pointer[j] = initialValue;
    } /* jForLoop */
    return(pointer) ;
}

// standard C memory allocation for our Grid Struct
double  *AllocateGridMemory ()
{
    int  j;
    double  *pointer;
    pointer = (double *)malloc(size * sizeof(Grid));
    for (j = 0; j < size; j++) {
        pointer[j] = initialValue;
    } /* jForLoop */
    return(pointer) ;
}

void  InitializeFdtd (Grid *g, int metalChoice)
{
    // Permittivity here is epsilon infinity now that everything is a Drude material
    double tempPermittivity, tempConductivity, tempPermeability, tempResistivity;
    double wPlasma, gDamping;

    // Switch block to determine what metal we're using:
    switch(metalChoice)
    {
      case 0: // Aluminum
        tempPermittivity = 1.0;
        tempConductivity = 1.0e+7;
        tempPermeability = 1.0;
        tempResistivity = 0.0;
        wPlasma = 1e+15; // WRONG
        gDamping = 0.1; // WRONG
        break;

      case 1: // Gold
        tempPermittivity = 1.0;
        tempConductivity = 1.0e+7;
        tempPermeability = 1.0;
        tempResistivity = 0.0;
        wPlasma = 1e+15; // WRONG
        gDamping = 0.1; // WRONG
        break;

      case 2: // Silver
        tempPermittivity = 1.0;
        tempConductivity = 1.0e+7;
        tempPermeability = 1.0;
        tempResistivity = 0.0;
        wPlasma = 1e+15; // WRONG
        gDamping = 0.1; // WRONG
        break;

      case 3: // Copper
        tempPermittivity = 1.0;
        tempConductivity = 1.0e+7;
        tempPermeability = 1.0;
        tempResistivity = 0.0;
        wPlasma = 1e+15; // WRONG
        gDamping = 0.1; // WRONG
        break;

    } /* Switch Block */

    /* Added for Drude metals so we can treat both the vacuum and object as "Drude"
       materials */
    double  mediaPlasma[MEDIACONSTANT] = {0.0, wPlasma}; // Plasma frequency
    double  mediaDamping[MEDIACONSTANT] = {0.0, gDamping}; // Damping constant
    /* End of Drude Metal Addition */

    double  mediaPermittivity[MEDIACONSTANT] = {1.0, tempPermittivity};    // eps, index=0 is for vacuum, index=1 is for the metallic cylinder
    double  mediaConductivity[MEDIACONSTANT] = {0.0, tempConductivity}; // sig,
    double  mediaPermeability[MEDIACONSTANT] = {1.0, tempPermeability};    // mur
    double  mediaResistivity[MEDIACONSTANT] = {0.0, tempResistivity};     // sim
    double  mediaCa[MEDIACONSTANT];
    double  mediaCb[MEDIACONSTANT];
    double  mediaDa[MEDIACONSTANT];
    double  mediaDb[MEDIACONSTANT];
    double  pi,speedOfLight,magneticPermeability0,electricalPermittivity0,frequency,wavelength,angularFrequency;
    double  dx,dt,reflectionCoefficient0,gradingOrder,temporary,electricalImpedance0,courantS,temp1,temp2,temp3;
    int  i,j,k, boundaryDataSize, media, boundaryIndex,xSizeMain,ySizeMain,numFreqs;
    int  abcSize ;
    double  cylinderDiameter, cylinderRadius, temporaryi,temporaryj,distance2 ;
    int  cylinderXCenter,cylinderYCenter;
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

    /***********************************************************************/
    //     Grid parameters
    /***********************************************************************/

    xSizeMain = 100;                              // number of main grid cells in x-direction
    ySizeMain = 50;                               // number of main grid cells in y-direction
    abcSize = ABCSIZECONSTANT;                    // thickness of PML region
    xSize = xSizeMain + 2 * abcSize;              // number of total grid cells in x-direction
    ySize = ySizeMain + 2 * abcSize;              // number of total grid cells in y-direction

    boundaryDataSize  = 2 * xSize * abcSize;                      // front edge + back edge
    boundaryDataSize += 2 * (abcSize * (ySize - 2 * abcSize));    // left + right edges

    //xSource = 50 + abcSize;                          //location of z-directed hard source
    //ySource = 50 + abcSize;                          //location of z-directed hard source
    xSource = 15 + abcSize;                       //location of z-directed hard source
    ySource = ySize / 2;                          //location of z-directed hard source

    dx = 1.0e-8;                                  //space increment of square lattice  (meters)
    dt = dx / (2.0 * speedOfLight);               //time step,  seconds, courant limit, Taflove1995 page 177
    courantS = dt / dx;                           // Courant Number

    maximumIteration = NUMBEROFITERATIONCONSTANT;                 //total number of time steps

    reflectionCoefficient0 = 1.0e-7;              // for PML, Nikolova part4 p.25
    gradingOrder = 3;                             // for PML, (m) was 2;  optimal values: 2 <= m <= 6,  Nikolova part4 p.29

    /***********************************************************************/
    //     Material parameters
    /***********************************************************************/

    media = MEDIACONSTANT;        // number of different medias, ie 2: vacuum, metallicCylinder

    /***********************************************************************/
    //     Wave excitation
    /***********************************************************************/

#if 0
    for (i = 0; i < maximumIteration; i++) {
        temporary = (double  )i;
        sourceValue[i] = sin( angularFrequency * temporary * dt);   // simple sine wave
    } /* iForLoop */
#endif
#if 1
    rtau = 160.0e-12;
    tau = rtau / dt;
    delay = 3 * tau;
    for (i = 0; i < maximumIteration; i++) {
        sourceValue[i] = 0.0;
    } /* iForLoop */
    for (i = 0; i < (int  )(7.0 * tau); i++) {
        temporary = (double  )i - delay;
        sourceValue[i] = sin( angularFrequency * (temporary) * dt) * exp(-( (temporary * temporary)/(tau * tau) ) );
    } /* iForLoop */
#endif

    /***********************************************************************/
    //     Field arrays
    /***********************************************************************/

    ex = AllocateMemory(xSize,    ySize + 1, 0.0 );        // 1 extra in y direction for pec
    ey = AllocateMemory(xSize + 1,ySize,     0.0 );        // 1 extra in x direction for pec
    hz = AllocateMemory(xSize,    ySize,     0.0 );
    hzy = AllocateMemory1D(boundaryDataSize, 0.0 );        // for the split-field data for hz in the pml regions

    e2Field = AllocateMemory(xSize, ySize, 0.0);           // E^2 for plotting

    /* Polarization Current Fields */
    jx = AllocateMemory(xsize, ySize, 0.0);
    jy = AllocateMemory(xsize, ySize, 0.0);
    exOld = AllocateMemory(xSize, ySize + 1, 0.0);
    eyOld = AllocateMemory(xSize + 1, ySize, 0.0);

    /* Pre-compute values that will be placed in to cjj and cje */
    double cjjTemp[media], cjeTemp[media];

    for (i = 0; i < media; i++) {
        cjjTemp[i] = (1 - (mediaDamping[i] * dt)/2) / (1 + (mediaDamping[i] * dt)/2);
        cjeTemp[i] = (1 / (1 + (mediaDamping[i] * dt)/2))) * (electricalPermittivity0*mediaPlasma[i]*mediaPlasma[i]*dt/2);
    } /* iForLoop */

    /***********************************************************************/
    //     DFT Array Initialization
    /***********************************************************************/
    numFreqs = NUMBERDFTFREQS;
    reflDFT = AllocateMemory1D(numFreqs, 0.0);
    tranDFT = AllocateMemory1D(numFreqs, 0.0);
    reflXPos = abcSize + 1;
    tranXPos = xSize - 1;

    double waveMin, waveMax; // Min / Max wavelength (expected in nm)
    waveMin = 400.0;
    waveMax = 800.0;

    // Find frequencies:
    for (i = 0; i < numFreqs; i++) {
      // Calculate evenly spaced frequencies from 400 to 800 nm
      temporary = waveMin/dx + ((double  )i * (waveMax - waveMin) / ((double  )numFreqs - 1));
      kList[i] = (int  )round(maximumIteration * courantS / temporary); // Schneider 5.29

      // Now convert back to find the "real" wavelength we are working with:
      wavelengthList[i] = maximumIteration * courantS / (double  )kList[i]; // Same, but inverted for wavelength

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
        temp3 = electricalImpedance0 * courantS / electricalPermittivity0;
        mediaCa[i] = (1.0 - temp1 - temp2) / (1.0 + temp1 + temp2);                                          // See Schneider 10.62
        mediaCb[i] = temp3 / (1.0 + temp1 + temp2);                                                          // ""

        temporary  = dt *  mediaResistivity[i] / (2.0 * magneticPermeability0 * mediaPermeability[i]);       // ditto
        mediaDa[i] = (1.0 - temporary) / (1.0 + temporary);                                                  // ditto
        mediaDb[i] = dt / (magneticPermeability0 * mediaPermeability[i] * dx * (1.0 + temporary));           // ditto
    } /* iForLoop */


    /***********************************************************************/
    //     Grid Coefficients
    /***********************************************************************/

    //     Initialize entire grid to free space

    caex = AllocateMemory(xSize, ySize, mediaCa[0]);     // note: don't need to allocate for pec region, as it is not evaluated
    cbex = AllocateMemory(xSize, ySize, mediaCb[0]);     // also: Initialize the entire grid to vacuum.
    caey = AllocateMemory(xSize, ySize, mediaCa[0] );
    cbey = AllocateMemory(xSize, ySize, mediaCb[0] );
    dahz = AllocateMemory(xSize, ySize, mediaDa[0] );
    dbhz = AllocateMemory(xSize, ySize, mediaDb[0] );
    dahzy = AllocateMemory1D(boundaryDataSize, mediaDa[0] );        // for the split-field data for hz in the pml regions
    dbhzy = AllocateMemory1D(boundaryDataSize, mediaDb[0] );        // for the split-field data for hz in the pml regions

    /* Initialize polarization current */
    cjj = AllocateMemory(xSize, ySize, 1.0);  // Initialize to 1 for free space
    cje = AllocateMemory(xSize, ySize, 0.0);  // Initialize to 0 for free space


    //     Add metal cylinder
#if 1
    cylinderDiameter = 20;                                     // diameter of cylinder
    cylinderRadius = cylinderDiameter / 2.0;                   // radius of cylinder
    cylinderXCenter = (4 * xSizeMain) / 5 + abcSize;           // i-coordinate of cylinder's center
    cylinderYCenter = ySize / 2;                               // j-coordinate of cylinder's center
    for (i = 0; i < xSize; i++) {
        for (j = 0; j < ySize; j++) {
            temporaryi = (double  )(i - cylinderXCenter) ;
            temporaryj = (double  )(j - cylinderYCenter) ;
            distance2 = (temporaryi + 0.5) * (temporaryi + 0.5) + (temporaryj) * (temporaryj);
            if (distance2 <= (cylinderRadius * cylinderRadius)) {
                caex[i][j] = mediaCa[1];
                cbex[i][j] = mediaCb[1];

                /* Polarization Current Constants: */
                cjj[i][j] = cjjTemp[1];
                cje[i][j] = cjeTemp[1];
            } /* if */
            // This looks tricky! Why can't caey/cbey use the same 'if' statement as caex/cbex above ??
            distance2 = (temporaryj + 0.5) * (temporaryj + 0.5) + (temporaryi) * (temporaryi);
            if (distance2 <= (cylinderRadius * cylinderRadius)) {
                caey[i][j] = mediaCa[1];
                cbey[i][j] = mediaCb[1];

                /* Polarization Current Constants: */
                cjj[i][j] = cjjTemp[1];
                cje[i][j] = cjeTemp[1];
            } /* if */
        } /* jForLoop */
    } /* iForLoop */
#endif

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


    /***********************************************************************/
    //     Print variables (diagnostic)
    /***********************************************************************/
/*
    printf("# pi:%16.10g\n",pi);
    printf("# speedOfLight:%16.10g\n",speedOfLight);
    printf("# magneticPermeability0:%16.10g\n",magneticPermeability0);
    printf("# electricalPermittivity0:%16.10g\n",electricalPermittivity0);
    printf("# frequency:%16.10g\n",frequency);
    printf("# wavelength:%16.10g\n",wavelength);
    printf("# angularFrequency:%16.10g\n",angularFrequency);
    printf("# dx:%16.10g\n",dx);
    printf("# dt:%16.10g\n",dt);
    printf("# reflectionCoefficient0:%16.10g\n",reflectionCoefficient0);
    printf("# gradingOrder:%16.10g\n",gradingOrder);
    printf("# xSizeMain:%d\n",xSizeMain);
    printf("# ySizeMain:%d\n",ySizeMain);
    printf("# abcSize:%d\n",abcSize);
    printf("# xSize:%d\n",xSize);
    printf("# ySize:%d\n",ySize);
    printf("# boundaryDataSize:%d\n",boundaryDataSize);
    printf("# xSource:%d\n",xSource);
    printf("# ySource:%d\n",ySource);
    printf("# maximumIteration:%d\n",maximumIteration);
*/

#if 0
    for (i = 0; i < maximumIteration; i++) {
        printf("%d: source: %16.10g\n",i,sourceValue[i]);
    } /* iForLoop */

    // print main grid geometry
    printf("total grid:\n");
    for (i = 0; i < xSize; i++) {
        printf("%3d: ",i);
        for (j = 0; j < ySize; j++) {
            ch = '.';
            if (caex[i][j] != mediaCa[0]) {
                ch = 'x';
            } /* if */
            if (caey[i][j] != mediaCa[0]) {
                if (ch == '.') {
                    ch = 'y';
                } /* if */
                else {
                    ch = 'B';
                } /* else */
            } /* if */
            if ((i == xSource) && (j == ySource)) {
                ch = 'S';
            } /* if */
            printf("%c",ch);
        } /* jForLoop */
        printf("\n");
    } /* iForLoop */
    printf("\n");

    for (j = ySize-1; j >= 0; j--) {
        printf("%3d: ",j);
        for (i = 0; i < xSize; i++) {
            printf("[%3d: ",i);
            // printf("ex:%+9.6g ",caex[i][j]);
            // printf("ey:%+9.6g ",caey[i][j]);
            // printf("hzx:%+9.6g ",dahz[i][j]);

            if ((i >= regionData[1].xStart) && (i < regionData[1].xStop) && (j >= regionData[1].yStart) && (j < regionData[1].yStop) ) {
                boundaryIndex = (j - regionData[1].yStart) + (i - regionData[1].xStart) * (regionData[1].yStop - regionData[1].yStart);
            } /* if */
            else if ((i >= regionData[2].xStart) && (i < regionData[2].xStop) && (j >= regionData[2].yStart) && (j < regionData[2].yStop) ) {
                boundaryIndex = (j - regionData[2].yStart) + (i - regionData[2].xStart) * (regionData[2].yStop - regionData[2].yStart) + (xSize * abcSize);
            } /* if */
            else {
                boundaryIndex = boundaryDataSize - 1;   // hack
            } /* else */
            printf("hzx:%+9.6g ",dahzy[boundaryIndex]);

            printf("] ");
        } /* jForLoop */
        printf("\n");
    } /* iForLoop */
    printf("\n");

    exit(0);
#endif

}




void  EvaluateFdtd (Grid *g, double  minimumValue, double  maximumValue)
{
    // Filename for where we will either save or load empty run data
    char filename[1000];
    FILE *filePointer;

    int  plottingInterval,iValue;        // for plotting
    int  n,i,j;
    int  boundaryIndex,regionIndex,xStart,yStart,xStop,yStop;
    double  hzx;                          // for Hz PML
    double  scaleValue,temporary;         // for plotting
    int  centerx = 50+ABCSIZECONSTANT ;   // for printing
    int  centery = 25+ABCSIZECONSTANT ;   //    ""

    /***********************************************************************/
    //     BEGIN TIME-STEPPING LOOP
    /***********************************************************************/

    plottingInterval = 0;
    for (n = 0; n  < maximumIteration; n++) {  // iteration loop

        fprintf(stderr,"n:%d\n",n);

        /***********************************************************************/
        //     Update electric fields (EX and EY)
        /***********************************************************************/

        for (i = 0; i < xSize; i++) {
            for (j = 1; j < ySize; j++) {        // j=0 = pec, so don't evaluate
                //ex[i][j] = caex[i][j] * ex[i][j] + cbex[i][j] * ( hz[i][j] - hz[i][j-1] ); /* Original update */
                exOld[i][j] = ex[i][j]; // Store previous field for polarization current
                temporary = (1.0/2.0) * (1.0 + cjj[i][j]) * jx[i][j];
                // If things go sideways, it might be the dx term on the H field... (same for y below)
                ex[i][j] = caex[i][j] * ex[i][j] + cbex[i][j] * ( dx * (hz[i][j] - hz[i][j-1]) - temporary ); // Need dx on H fields?
            } /* jForLoop */
        } /* iForLoop */

        for (i = 1; i < xSize; i++) {            // i=0 = pec, so don't evaluate
            for (j = 0; j < ySize; j++) {
                //ey[i][j] = caey[i][j] * ey[i][j] + cbey[i][j] * ( hz[i-1][j] - hz[i][j] ); /* Original update */
                eyOld[i][j] = ey[i][j]; // Store previous field for polarization current
                temporary = (1.0/2.0) * (1.0 + cjj[i][j]) * jy[i][j];
                ey[i][j] = caey[i][j] * ey[i][j] + cbey[i][j] * ( dx * (hz[i-1][j] - hz[i][j]) - temporary ); // Need dx on H fields?
            } /* jForLoop */
        } /* iForLoop */


        /***********************************************************************/
        //     Update magnetic fields (HZ) in center (main) grid
        /***********************************************************************/


        regionIndex = 0;    // center (main) grid
        xStart = regionData[regionIndex].xStart;
        xStop  = regionData[regionIndex].xStop ;
        yStart = regionData[regionIndex].yStart;
        yStop  = regionData[regionIndex].yStop ;
        for (i = xStart; i < xStop; i++) {
            for (j = yStart; j < yStop; j++) {
                hz[i][j] = dahz[i][j] * hz[i][j] + dbhz[i][j] * ( ex[i][j+1] - ex[i][j] + ey[i][j] - ey[i+1][j] );
            } /* jForLoop */
        } /* iForLoop */

        hz[xSource][ySource] = sourceValue[n];


        /***********************************************************************/
        //     Update HZ in PML regions (hzx,hzy)
        /***********************************************************************/

        boundaryIndex = 0;
        for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
            xStart = regionData[regionIndex].xStart;
            xStop  = regionData[regionIndex].xStop ;
            yStart = regionData[regionIndex].yStart;
            yStop  = regionData[regionIndex].yStop ;
            for (i = xStart; i < xStop; i++) {
                for (j = yStart; j < yStop; j++) {
                    hzx = hz[i][j] - hzy[boundaryIndex];   // extract hzx
                    hzx = dahz[i][j] * hzx + dbhz[i][j] * ( ey[i][j] - ey[i+1][j] );    // dahz,dbhz holds dahzx,dbhzx
                    hzy[boundaryIndex] = dahzy[boundaryIndex] * hzy[boundaryIndex] + dbhzy[boundaryIndex] * ( ex[i][j+1] - ex[i][j] );
                    hz[i][j] = hzx +  hzy[boundaryIndex];  // update hz
                    boundaryIndex++;
                } /* jForLoop */
            } /* iForLoop */
        } /* forLoop */


        /***********************************************************************/
        //     Update polarization current (jx,jy) (actually are dx*jx or dx*jy)
        /***********************************************************************/
        for (i = 0; i < xSize; i++) {
            for (j = 1; j < ySize; j++) {        // j=0 = pec, so don't evaluate
                jx[i][j] = cjj[i][j] * dx * jx[i][j] + cje[i][j] * (ex[i][j] + exOld[i][j]);
            } /* jForLoop */
        } /* iForLoop */

        for (i = 1; i < xSize; i++) {            // i=0 = pec, so don't evaluate
            for (j = 0; j < ySize; j++) {
                jy[i][j] = cjj[i][j] * dx * jy[i][j] + cje[i][j] * (ey[i][j] + eyOld[i][j]);
            } /* jForLoop */
        } /* iForLoop */


        /***********************************************************************/
        //     Update DFT values
        /***********************************************************************/
        regionIndex = 0;    // center (main) grid
        xStart = regionData[regionIndex].xStart;
        xStop  = regionData[regionIndex].xStop ;
        yStart = regionData[regionIndex].yStart;
        yStop  = regionData[regionIndex].yStop ;

        for (i = 0; i < NUMBERDFTFREQS; i++) {
          for (j = yStart; j < yStop; j++) {
            temporary = -1.0 * 0.5 * (ey[reflXPos][j]*hz[reflXPos][j]); // Poynting Flux through our line (positive x) (so we negate the result)
            reflDFT[i] += temporary * cos(2*pi*kList[i]*n/maximumIteration); // Schneider 5.32

            temporary = 0.5 * (ey[tranXPos][j]*hz[tranXPos][j); // Poynting Flux through our line
            tranDFT[i] += temporary * cos(2*pi*kList[i]*n/maximumIteration);
          } /* jForLoop */
        } /* iForLoop */

        /***********************************************************************/
        //     Plot fields
        /***********************************************************************/
#if 1
        if (plottingInterval == 0) {
            plottingInterval = 2;

            // Send E^2 field to our plotting routine:
            scaleValue = 256.0 / (maximumValue - minimumValue);
            for (j = 0; j < ySize; j++) {
                for (i = 0; i < xSize; i++) {
                    e2Field[i][j] = ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j];
                } /* xForLoop */
            } /* yForLoop */
            imageShow(e2Field);

        } /* if */
        plottingInterval--;
#endif

    } /* iteration forLoop */

// Switch for empty run to export data vs real run
#if 1 // Empty run

    // Write to header file for use later
    reflFile = "empty_refl_data.h";
    reflPointer = fopen(reflFile, "w");

    tranFile = "empty_tran_data.h";
    tranPointer = fopen(filename, "w");

    for (i = 0; i < NUMBERDFTFREQS; i++) {
      reflDFT[i] = reflDFT[i] / maximumIterations;
      tranDFT[i] = tranDFT[i] / maximumIterations;

      putc(reflDFT[i], reflPointer);
      putc(tranDFT[i], tranPointer);
    } /* iForLoop */

    fclose(reflPointer);
    fclose(tranPointer);

#endif

#if 0 // Real run

    // Normalize our DFT's by number of steps taken and relative to empty run:
    for (i = 0; i < NUMBERDFTFREQS; i++) {
      reflDFT[i] = reflDFT[i] / ( maximumIterations * emptyReflDFT[i] );
      tranDFT[i] = tranDFT[i] / ( maximumIterations * emptyTranDFT[i] );
    } /* iForLoop */

#endif

}


int  main ()
{
    // Allocate Grid Memory:
    Grid *g;
    g = AllocateGridMemory();

    InitializeFdtd(g, 0);
    EvaluateFdtd(g, -0.1, 0.1);
    return(0) ;
}


