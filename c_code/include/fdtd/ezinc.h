#ifndef EZINC
#define EZINC

#include "fdtd_macro.h"


double ezInc(double time, double location);
void ezIncInit(struct Grid *g);
void lineSource(struct Grid *g, int x_ind, int time);

#endif
