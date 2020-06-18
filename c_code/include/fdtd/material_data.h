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
#define CONVERSION 6.28319/PLANCK_H // Conversion for constants that have units of energy (2 Pi/h) to get to rad/sec

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
    .epsInf = 1.0,
    .permeability = 1.0,
    .conductivity = 413.505 * EPS0 * CONVERSION,
    .params[0] = {
      .ap = (-0.272493 + I*1.44484)*CONVERSION,
      .cp = (3.80577 - I*5.88315)*CONVERSION,
    },
    .params[1] = {
      .ap = (-0.116623 + I*0.133088)*CONVERSION,
      .cp = (-209.786 - I*295.234)*CONVERSION,
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
    .epsInf = 3.35141,
    .permeability = 1.0,
    .conductivity = 105.955 * EPS0 * CONVERSION,
    .params[0] = {
      .ap = (-0.063629 + I*0.154755)*CONVERSION,
      .cp = (-55.7092 - I*260.36)*CONVERSION,
    },
    .params[1] = {
      .ap = (5.26744 - I*3.11234)*CONVERSION,
      .cp = (-17.5676 - I*34.0369)*CONVERSION,
    },
  },
  { //3, Cu, J&C
    .num_poles = 2,
    .epsInf = 1.98801,
    .permeability = 1.0,
    .conductivity = 27.8466 * EPS0 * CONVERSION,
    .params[0] = {
      .ap = (-0.607993 + I*2.0456)*CONVERSION,
      .cp = (3.84516 + I*1.37118)*CONVERSION,
    },
    .params[1] = {
      .ap = (-0.0173655 - I*0.134845)*CONVERSION,
      .cp = (-9.84239 + I*261.159)*CONVERSION,
    },
  },
  { //4, SiO2
    .num_poles = 0,
    .epsInf = 2.2201, // constant n = 1.49
    .permeability = 1.0,
    .conductivity = 0.0,
  },
  { //5, Si, Green, 2008
    .num_poles = 3,
    .epsInf = 1.30676,
    .permeability = 1.0,
    .conductivity = 0.0471157 * EPS0 * CONVERSION,
    .params[0] = {
      .ap = (-0.116424 + I*3.3587)*CONVERSION,
      .cp = (1.96042 - I*2.20821)*CONVERSION,
    },
    .params[1] = {
      .ap = (-0.381842 + I*4.26879)*CONVERSION,
      .cp = (-0.578006 - I*13.8223)*CONVERSION,
    },
    .params[2] = {
      .ap = (-0.37732 - I*3.6481)*CONVERSION,
      .cp = (1.38221 + I*4.59924)*CONVERSION,
    },
  },
};

#endif
