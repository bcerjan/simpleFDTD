#include "ezinc.h"
#include <math.h>

static double cdtds, ppw = 0;

/* Initialize Variables for source */
void ezIncInit(struct Grid *g) {
  ppw = 6e-7 / dx; // Fixed source at 600 nm wavlength
  cdtds = courantS;

  return;
}

/* Calculate source function at given time and location */
double ezInc(double time, double location) {
  double arg = M_PI*((cdtds * time - location) / ppw - 1.0); // magic 800 is bad, figure out why this is wrong...
  arg = arg * arg;

  return (1.0 - 2.0*arg)*exp(-arg);
}

/* Function to set one edge of the simulation (in x) to be equal to source term */
void lineSource(struct Grid *g, int x_ind, int time) {
  int j;
  double t = time;
  //if(time < 276) {
    for (j = 1; j < ySize-1; j++) {
      ey[x_ind][j] += ezInc(t, 0.0); // Y-Polarized source field
    } /* jForLoop */
  //} /* if */
  return;
}
