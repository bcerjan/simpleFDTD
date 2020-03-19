#ifndef SDL_IMAGE_FUNC
#define SDL_IMAGE_FUNC

#include "fdtd_macro.h"

double **findMatEdge(double **matrix, int sizeX, int sizeY);
void imageShow(struct Grid *g);
void PlotField (struct Grid *g, double maximumValue, double minimumValue);

#endif
