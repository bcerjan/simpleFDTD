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
#include "fdtd_grid.h"

// standard C memory allocation for 2-D array
double **AllocateMemory (int imax, int jmax, double initialValue)
{
    int  i,j;
    double  **pointer;
    pointer = malloc(imax * sizeof(double *));
    if (pointer == NULL) {
        printf("Error! memory not allocated.\n");
        exit(0);
    }

    for (i = 0; i < imax; i++) {
        pointer[i] = malloc(jmax * sizeof(double));
        if (pointer[i] == NULL) {
            printf("Error! memory not allocated.\n");
            exit(0);
        }
        for (j = 0; j < jmax; j++) {
            pointer[i][j] = initialValue;
        } /* jForLoop */
    } /* iForLoop */
    return(pointer) ;
}

// standard C memory allocation for 1-D array
double *AllocateMemory1D (int size, double initialValue)
{
//    printf("In 1D Allocation...\n");
    int  j;
    double  *pointer;
    pointer = malloc(size * sizeof(double));
    if (pointer == NULL) {
        printf("Error! memory not allocated.\n");
        exit(0);
    }

//    printf("1D pointer generated...\n");
    for (j = 0; j < size; j++) {
        pointer[j] = initialValue;
    } /* jForLoop */
    return(pointer) ;
}

// standard memory allocation for 3-D array
double ***AllocateMemory3D (int  imax, int  jmax, int kmax, double *initArray) {
  int i;
  double ***pointer;
  pointer = malloc(imax * sizeof(double **));
  if (pointer == NULL) {
      printf("Error! memory not allocated.\n");
      exit(0);
  }
  for (i = 0; i < imax; i++) {
    pointer[i] = AllocateMemory(jmax, kmax, initArray[i]);
  } /* iForLoop */
  return(pointer);
}

// standard C memory allocation for our Grid Struct
struct Grid *AllocateGridMemory ()
{
    struct Grid  *pointer = malloc(sizeof(struct Grid));

    if (pointer == NULL) {
        printf("Error! memory not allocated.\n");
        exit(0);
    }
    return(pointer) ;
}
