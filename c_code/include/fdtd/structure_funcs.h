#ifndef STRUCT_FUNCS
#define STRUCT_FUNCS

#include "fdtd_macro.h"

void structInit(int xCenter, int yCenter);
void addDisk(struct Grid *g, double radius);
void addRect(struct Grid *g, double width, double length);
void addTriangle(struct Grid *g, double length);

#endif
