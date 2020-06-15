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


/** Static definition of our material data to clean up grid initialization and
    make it easier to add new materials

    This data was taken from the source (typically via refractiveindex.info),
    truncated to the range we are getting results for (200-1000 nm), and
    then fit with a critical-points model (simultaneously for Re and Im components):
      e[w] = e_Inf + sigma / (-i * e_0 * w) + Sum[ -c_p / (i*w + a_p) + c.c ]

      (see Prokopidis and Zografopoulos eq. 1)

   Current list of materials: 0 -> Al, 1 -> Au, 2 -> Ag, 3 -> Cu, 4 -> Silica,
   5 -> Silicon,
**/


#ifndef MATERIAL_DATA
#define MATERIAL_DATA

// Maximumum number of poles we'll ever potentially need
#define MAX_POLES (8)

#include <complex.h>

#define PLANCK_H 4.135667e-15 // eV/sec
#define EPS0 8.854187e-12
#define SPEEDC 2.99792458e8 // m/sec
#define CONVERSION 6.28319/PLANCK_H // Conversion for constants that have units of energy (2 PI/h) to get to rad/sec

struct cpParams {
  complex double ap;
  complex double cp;
};

struct Material {
  int num_poles;
  double epsInf;
  double permeability;
  double conductivity;

  struct cpParams params[MAX_POLES];
};



static const struct Material materialData[6] = {
  { //0, Al, Vial, 2011
    .num_poles = 2,
    .epsInf = 0.329666,
    //.epsInf = 1.329666,
    .permeability = 1.0,
    .conductivity = 254.536 * EPS0 * CONVERSION,
    .params[0] = {
      .ap = (-0.0783386 + I*0.0699538)*CONVERSION,
      .cp = (-128.853 - I*732.243)*CONVERSION,
    },
    .params[1] = {
      .ap = (-0.255539 - I*1.46655)*CONVERSION,
      .cp = (2.81621 + I*5.36337)*CONVERSION,
    },
  },
  { //1, Au, Johnson and Christy
    .num_poles = 2,
    .epsInf = 1.0,
    .permeability = 1.0,
    .conductivity = 5.19039 * EPS0 * CONVERSION,
    .params[0] = {
      .ap = (-0.994406 + I*2.45951)*CONVERSION,
      .cp = (7.01819 - I*0.31594)*CONVERSION,
    },
    .params[1] = {
      .ap = (-0.0219569 + I*0.224007)*CONVERSION,
      .cp = (-0.67478 - I*155.767)*CONVERSION,
    },
  },
  { //2, Ag, Wu, 2014
    .num_poles = 2,
    .epsInf = 2.17741,
    .permeability = 1.0,
    .conductivity = 1.0e12,
    .params[0] = {
      .ap = (-0.388929 - I*2.22844)*CONVERSION,
      .cp = (4.26389 + I*8.17685)*CONVERSION,
    },
    .params[1] = {
      .ap = (-0.388929 - I*2.22844)*CONVERSION,
      .cp = (4.26389 + I*8.17685)*CONVERSION,
    },
  },
  { //3, Cu, J&C
    .num_poles = 2,
    .epsInf = 1.82316,
    .permeability = 1.0,
    .conductivity = 1.0e12,
    .params[0] = {
      .ap = (-0.388929 - I*2.22844)*CONVERSION,
      .cp = (4.26389 + I*8.17685)*CONVERSION,
    },
    .params[1] = {
      .ap = (-0.388929 - I*2.22844)*CONVERSION,
      .cp = (4.26389 + I*8.17685)*CONVERSION,
    },
  },
  { //4, SiO2, Lemarchand, 2013
    .num_poles = 0,
    .epsInf = 2.22439,
    .permeability = 1.0,
    .conductivity = 1.0e-10, // needs to be non-0 for stability
  },
  { //5, Si, Green, 2008
    .num_poles = 3,
    .epsInf = 0.793124,
    .permeability = 1.0,
    .conductivity = 1.0e12,
    .params[0] = {
      .ap = (-0.388929 - I*2.22844)*CONVERSION,
      .cp = (4.26389 + I*8.17685)*CONVERSION,
    },
    .params[1] = {
      .ap = (-0.388929 - I*2.22844)*CONVERSION,
      .cp = (4.26389 + I*8.17685)*CONVERSION,
    },
    .params[2] = {
      .ap = (-0.388929 - I*2.22844)*CONVERSION,
      .cp = (4.26389 + I*8.17685)*CONVERSION,
    },
  },
};

#endif
