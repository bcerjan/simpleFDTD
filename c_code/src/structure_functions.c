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

/* File containing functions to define different geometries for structure */
#include <math.h>
#include "structure_funcs.h"

static int xCent,yCent;

void structInit(int xCenter, int yCenter){
  xCent = xCenter;
  yCent = yCenter;
  return;
}

/* Add a disk. Radius is given in meters */
void addDisk(struct Grid *g, double xRadius, double yRadius) {
  int i,j;
  double tempi, tempj, distanceX, distanceY;
  for (i = 0; i < xSize; i++) {
    for (j = 0; j < ySize; j++) {
      tempi = (double )(i - xCent) * dx;
      tempj = (double )(j - yCent) * dx;
      distanceX = (tempi*tempi)/(xRadius*xRadius); // This is offset from "perfect" by half a unit cell
      distanceY = (tempj*tempj)/(yRadius*yRadius);
      if(distanceX + distanceY <= 1.0) {
        object_locs[i][j] = 1.0;
      } /* if */
    } /* jForLoop */
  } /* iForLoop */

  return;
}

/* Add a rectangle. Width is in x direction, length is in y. Dimensions in meters */
void addRect(struct Grid *g, double width, double length) {
  int i,j;
  double tempi, tempj, distanceX,distanceY;

  for (i = 0; i < xSize; i++) {
    for (j = 0; j < ySize; j++) {
      tempi = (double )(i - xCent) * dx;
      tempj = (double )(j - yCent) * dx;
      distanceX = fabs(tempi); // This is offset from "perfect" by half a unit cell
      distanceY = fabs(tempj); // ""
      if(distanceX <= width/2.0 && distanceY <= length/2.0) {
        object_locs[i][j] = 1.0;
      } /* if */
    } /* jForLoop */
  } /* iForLoop */

  return;
}


/* Add a triangle. Length = total 'size' in a given direction. Pointed towards x = 0 */
void addTriangle(struct Grid *g, double xLength, double yLength) {
  int i,j;
  double tempi, tempj, line1, line2;

  //d1 = length / (2.0 * sin(30.0 * M_PI/180.0)); // Offset distance from center point in positive-x direction
  double pt1[2] = { (xCent * dx) - xLength/2.0, (yCent * dx) }; // Left point
  double pt2[2] = { (xCent * dx) + xLength/2.0, (yCent * dx) + yLength/2.0 }; // Top Point
  double pt3[2] = { (xCent * dx) + xLength/2.0, (yCent * dx) - yLength/2.0 }; // Bottom Point

  for (i = 0; i < xSize; i++) {
    for (j = 0; j < ySize; j++) {
      tempi = (double )(i) * dx;
      tempj = (double )(j) * dx;
      // Lines definining top (line1) and bottom (line2) of our triangular region
      line1 = pt1[1] + ( (pt2[1]-pt1[1]) / (pt2[0] - pt1[0]) ) * (tempi - pt1[0]);
      line2 = pt3[1] + ( (pt1[1]-pt3[1]) / (pt1[0] - pt3[0]) ) * (tempi - pt3[0]);

      if(tempi <= pt3[0] && tempj <= line1 && tempj >= line2 ) {
        object_locs[i][j] = 1.0;
      } /* if */
    } /* jForLoop */
  } /* iForLoop */


  return;
}
