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
