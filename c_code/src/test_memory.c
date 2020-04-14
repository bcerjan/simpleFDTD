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
#include "array_proto.h"
#include "fdtd_proto.h"
#include "fdtd_macro.h"

/** Function to find walls using arma mat instead of image: **/
void findMatEdge(struct Grid *g) {
  int i,j,n,k;

  //double **out = AllocateMemory(sizeX, sizeY, 0.0);
  double gxVal;
  double gyVal;
  double test = 0.0;
  double gx[3][3] = {0};
  double gy[3][3] = {0};

  // Specify Sobel Operators:
  gx[0][0] = -1.0;
  gx[0][2] = 1.0;
  gx[1][0] = -2.0;
  gx[1][2] = 2.0;
  gx[2][0] = -1.0;
  gx[2][2] = 1.0;

  gy[0][0] = 1.0;
  gy[0][1] = 2.0;
  gy[0][2] = 1.0;
  gy[2][0] = -1.0;
  gy[2][1] = -2.0;
  gy[2][2] = -1.0;

/*
  //Sobel Operators:
  mat gx = { {-1, 0, 1},
             {-2, 0, 2},
             {-1, 0, 1} };

  mat gy = { {1, 2, 1},
             {0, 0, 0},
             {-1, -2, -1} };
*/
  for (i = 1; i < xSize - 3; i++) {
    for (j = 1; j < ySize - 3; j++) {
      gxVal = 0.0;
      gyVal = 0.0;
      // Manual 2-D convolution:
      for (n = -1; n < 2; n++) { // CHECK ME. MIGHT BE OFF-BY-ONE!!
        for (k = -1; k < 2; k++) {
          gxVal += gx[n+1][k+1] * object_locs[i+n][j+k];
          gyVal += gy[n+1][k+1] * object_locs[i+n][j+k];
          if (j == ySize/2 && i == xSize - ABCSIZECONSTANT - 30) {
            printf("object_locs[%i][%i]: %f\n", i+n,j+k, object_locs[i+n][j+k] );
            printf("gxVal[%i][%i]: %f\n", n,k, gxVal );
            printf("gyVal[%i][%i]: %f\n", n,k, gyVal );
          }
        } /* kForLoop */
      } /* nForLoop */

      edgeMat[i+1][j+1] = pow(gxVal*gxVal + gyVal*gyVal, 0.5);
      /*if (j == ySize/2) {
        printf("edgeMat[%i][%i]: %f\n",i+1,j+1,edgeMat[i+1][j+2]);
      }*/
    } /* jForLoop */
  } /* iForLoop */

  return;

}



int main()
{

  struct Grid *g;
  g = AllocateGridMemory();

  InitializeFdtd(g, 0, 0, 100.0, 1.0, 1.0);

  findMatEdge(g);


  printf("Max edgeMat: %f\n", ArrayMax(edgeMat,xSize,ySize));
  printf("Min edgeMat: %f\n", ArrayMin(edgeMat,xSize,ySize));
  freeGrid(g);

  return(0);
}
