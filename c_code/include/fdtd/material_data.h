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
    truncated to the range we are getting results for (400-800 nm), and
    then fit with the critical-points model (simultaneously for Re and Im components):
      χ[f] = 	Sum[ Ap * Omegap * (exp[-i*phip]/(Omegap + f - i*Gammap) + exp[i*phip]/(Omegap - f + i*Gammap))

      (see Prokopidis and Dimitrios eq. 1 with no Drude term {it has been folded in
      to the other critical point terms, which wastes some memory, but simplifies
      the coding})

    Note that some of the fits are actually quite poor! This is just a rough
    idea of what a particular metal will do and is not a great substitute for
    a more elaborate simulation.

    Order is: epsilon infinity, conductivity, permeability, resistivity,
    plasma frequency, damping rate: (eps, sigma, mu, sim, f_p, g)
    values in meters, Hz - type units (not radians)

   Current list of materials: 0 -> Al, 1 -> Au, 2 -> Ag, 3 -> Cu, 4 -> Silica
**/


#ifndef MATERIAL_DATA
#define MATERIAL_DATA

// Maximumum number of poles we'll eveer potentially need
#define MAX_POLES (12)

struct cpParams {
  double bigA;
  double Omega;
  double phi;
  double Gamma;
};

struct Material {
  int num_poles;
  double epsInf;
  double permeability;
  double conductivity;
  struct cpParams params[MAX_POLES];
};



static const struct Material materialData[5] = {
  { //0, Al, Cheng 2016
    .num_poles = 2,
    .epsInf = 4.1,
    .permeability = 1.0,
    .conductivity = 1.0e12,
    .params[0] = {
      .bigA = 2.0,
      .Omega = 1.0,
      .phi = 0.01,
      .Gamma = -0.5,
    },
    .params[1] = {
      .bigA = 20.0,
      .Omega = 10.0,
      .phi = 00.01,
      .Gamma = -0.05,
    }
  }
};

#endif
