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
#include "fdtd_field_proto.h"

/* Function to run our FDTD simulation */

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
                ex[i][j] = caex[i][j] * ex[i][j] + cbex[i][j] * (dx * (hz[i][j] - hz[i][j-1]) - temporary ); // Need dx on H fields?
            } /* jForLoop */
        } /* iForLoop */

        for (i = 1; i < xSize; i++) {            // i=0 = pec, so don't evaluate
            for (j = 0; j < ySize; j++) {
                //ey[i][j] = caey[i][j] * ey[i][j] + cbey[i][j] * ( hz[i-1][j] - hz[i][j] ); /* Original update */
                eyOld[i][j] = ey[i][j]; // Store previous field for polarization current
                temporary = (1.0/2.0) * (1.0 + cjj[i][j]) * jy[i][j];
                ey[i][j] = caey[i][j] * ey[i][j] + cbey[i][j] * (dx * (hz[i-1][j] - hz[i][j]) - temporary ); // Need dx on H fields?
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
