#ifndef SDL_IMAGE_FUNC
#define SDL_IMAGE_FUNC

#include "fdtd_macro.h"

void imageInit(struct Grid *g);
void imageFree();
void findMatEdge(struct Grid *g);
void imageShow(struct Grid *g);
void PlotField (struct Grid *g, double maximumValue, double minimumValue);

#endif
