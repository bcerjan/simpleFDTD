#include <stdlib.h>
#include <stdio.h>
#include "fdtd_grid.h"

// standard C memory allocation for 2-D array
double  **AllocateMemory (int  imax, int  jmax, double  initialValue)
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
double  *AllocateMemory1D (int  size, double  initialValue)
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

// standard C memory allocation for our Grid Struct
struct Grid  *AllocateGridMemory ()
{
    struct Grid  *pointer = malloc(sizeof(struct Grid));

    if (pointer == NULL) {
        printf("Error! memory not allocated.\n");
        exit(0);
    }
    return(pointer) ;
}