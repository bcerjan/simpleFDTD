#ifndef FDTD_PROTOTYPES
#define FDTD_PROTOTYPES

#include "fdtd_grid.h"

void InitializeFdtd(struct Grid *g, int metalChoice, int objectChoice);
void freeGrid(struct Grid *g);
double **AllocateMemory(int imax, int jmax, double initialValue);
double *AllocateMemory1D(int size, double initialValue);
double *AllocateGridMemory();
void EvaluateFdtd(struct Grid *g, double minimumValue, double maximumValue);

#endif
