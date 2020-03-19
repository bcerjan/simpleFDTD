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
void addDisk(struct Grid *g, double radius) {
  int i,j;
  double tempi, tempj, distance2;
  for (i = 0; i < xSize; i++) {
    for (j = 0; j < ySize; j++) {
      tempi = (double )(i - xCent) * dx;
      tempj = (double )(j - yCent) * dx;
      distance2 = (tempi*tempi) + (tempj*tempj); // This is offset from "perfect" by half a unit cell
      if(distance2 <= (radius * radius)) {
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
      distanceX = fabs(tempi - width/2.0); // This is offset from "perfect" by half a unit cell
      distanceY = fabs(tempj - length/2.0); // ""
      if(distanceX <= width/2.0 && distanceY <= length/2.0) {
        object_locs[i][j] = 1.0;
      } /* if */
    } /* jForLoop */
  } /* iForLoop */

  return;
}


/* Add an equilateral triangle. Length = side length. Pointed towards x = 0 */
void addTriangle(struct Grid *g, double length) {
  int i,j;
  double tempi, tempj, line1, line2, d1, d2;

  d1 = length / (2.0 * sin(30.0 * M_PI/180.0)); // Offset distance from center point in positive-x direction
  d2 = length * sin(60.0 * M_PI/180.0) - d1;
  double pt1[2] = { (xCent * dx) + d1, (yCent * dx) + length/2.0 };
  double pt2[2] = { (xCent * dx) - d2, (yCent * dx) };
  double pt3[2] = { (xCent * dx) + d1, (yCent * dx) - length/2.0 };

  for (i = 0; i < xSize; i++) {
    for (j = 0; j < ySize; j++) {
      tempi = (double )(i) * dx;
      tempj = (double )(j) * dx;
      // Lines definining top (line1) and bottom (line2) of our triangular region
      line1 = pt1[1] + ( (pt2[1]-pt1[1]) / (pt2[0] - pt1[0]) ) * (tempi - pt1[0]);
      line2 = pt3[1] + ( (pt2[1]-pt3[1]) / (pt2[0] - pt3[0]) ) * (tempi - pt3[0]);

      if(tempi <= pt3[0] && tempj <= line1 && tempj >= line2 ) {
        object_locs[i][j] = 1.0;
      } /* if */
    } /* jForLoop */
  } /* iForLoop */


  return;
}
