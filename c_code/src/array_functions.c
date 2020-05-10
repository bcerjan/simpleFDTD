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

#include <math.h>
#include <stdio.h>
#include "fdtd_macro.h"
/* File containing array functions (max, min) */

double ArrayMax(double **ptr, int xWidth, int yWidth)
{
  int i,j;
  double max = ptr[0][0];
  for (i = 0; i < xWidth; i++) {
    for (j = 0; j < yWidth; j++) {
      if (ptr[i][j] > max) {
        max = ptr[i][j];
      } /* ifBlock */
    } /* jForLoop */
  } /* iForLoop */
  return max;
}

double ArrayMin(double **ptr, int xWidth, int yWidth)
{
  int i,j;
  double min = ptr[0][0];
  for (i = 0; i < xWidth; i++) {
    for (j = 0; j < yWidth; j++) {
      if (ptr[i][j] < min) {
        min = ptr[i][j];
      } /* ifBlock */
    } /* jForLoop */
  } /* iForLoop */
  return min;
}

double AbsArrayMax(double **ptr, int xWidth, int yWidth)
{
  int i,j,itemp=0,jtemp=0;
  double max = fabs(ptr[0][0]);
  double temp = 0.0;
  for (i = 0; i < xWidth; i++) {
    for (j = 0; j < yWidth; j++) {
      temp = fabs(ptr[i][j]);
      if (temp > max) {
        max = temp;
        itemp = i;
        jtemp = j;
      } /* ifBlock */
    } /* jForLoop */
  } /* iForLoop */
  printf("[i][j]: [%i][%i]\n",itemp,jtemp );
  return max;
}

double AbsArrayMin(double **ptr, int xWidth, int yWidth)
{
  int i,j;
  double min = fabs(ptr[0][0]);
  double temp = 0.0;
  for (i = 0; i < xWidth; i++) {
    for (j = 0; j < yWidth; j++) {
      temp = fabs(ptr[i][j]);
      if (temp < min) {
        min = temp;
      } /* ifBlock */
    } /* jForLoop */
  } /* iForLoop */
  return min;
}

double AbsVectorMax(double *ptr, int length) {
  int i,itemp;
  itemp = 0;
  double max = fabs(ptr[0]);
  double temp = 0.0;
  for (i = 0; i < length; i++) {
    temp = fabs(ptr[i]);
    if (temp > max) {
      max = temp;
      itemp = i;
    }
  }
  printf("[i]: [%i]\n",itemp);
  return max;
}

/* Function to find max value in PML region of a 1-D pointer defined there */
double AbsPML1DMax(struct Grid *g, double *ptr)
{
  int i,j,boundaryIndex,regionIndex,xStop,xStart,yStop,yStart;
  double max = 0.0;
  double temp = 0.0;
  boundaryIndex = 0;
  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;

    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
        temp = fabs(ptr[boundaryIndex]);
        if (temp > max) {
          max = temp;
        } /* ifBlock */
        boundaryIndex++;
      } /* jForLoop */
    } /* iForLoop */
  } /* region forLoop */
  return max;
}

/* Function to find max value in PML region of a 2-D pointer defined there */
double AbsPML2DMax(struct Grid *g, double **ptr)
{
  int i,j,regionIndex,xStop,xStart,yStop,yStart;
  double max = 0.0;
  double temp = 0.0;

  for (regionIndex = 1; regionIndex < NUMBEROFREGIONS; regionIndex++) {
    xStart = regionData[regionIndex].xStart;
    xStop  = regionData[regionIndex].xStop ;
    yStart = regionData[regionIndex].yStart;
    yStop  = regionData[regionIndex].yStop ;

    for (i = xStart; i < xStop; i++) {
      for (j = yStart; j < yStop; j++) {
        temp = fabs(ptr[i][j]);
        if (temp > max) {
          max = temp;
        } /* ifBlock */
      } /* jForLoop */
    } /* iForLoop */
  } /* region forLoop */
  return max;
}
