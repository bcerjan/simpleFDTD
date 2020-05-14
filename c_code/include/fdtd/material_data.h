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

// Maximumum number of poles we'll ever potentially need
#define MAX_POLES (8)

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
  double drudePlasma;
  double drudeDamping;
  struct cpParams params[MAX_POLES];
};



static const struct Material materialData[5] = {
  { //0, Al, Cheng 2016
    .num_poles = 3,
    .epsInf = 504.939,
    .permeability = 1.0,
    .conductivity = 1.0e12,
    .drudePlasma = 20.064e14,
    .drudeDamping = 2.48281e14,
    .params[0] = {
      .bigA = 223.167,
      .Omega = 27.1198e14,
      .phi = -6.15405,
      .Gamma = 2.11272e9,
    },
    .params[1] = {
      .bigA = 470.935,
      .Omega = 35.9801e14,
      .phi = -3.04908,
      .Gamma = 3.87203e9,
    },
    .params[2] = {
      .bigA = 1.84289,
      .Omega = 3.45883e14,
      .phi = -1.19693,
      .Gamma = 0.475579e14,
    }
  },
  { //2, Au, Calculated
    .num_poles = 2,
    .epsInf = 1.11781,
    .permeability = 1.0,
    .conductivity = 1.0e12,
    .drudePlasma = 13.1839e15,
    .drudeDamping = 0.109093e15,
    .params[0] = {
      .bigA = 3.0392,
      .Omega = 4.20979e15,
      .phi = -1.09027,
      .Gamma = 2.35485e15,
    },
    .params[1] = {
      .bigA = 0.273624,
      .Omega = 3.88228e15,
      .phi = -0.412269,
      .Gamma = 0.452509e15,
    },
  },
  { //2, Ag, Calculated
    .num_poles = 2,
    .epsInf = 2.17741,
    .permeability = 1.0,
    .conductivity = 1.0e12,
    .drudePlasma = 14.4203e15,
    .drudeDamping = 0.0262676e15,
    .params[0] = {
      .bigA = 0.855074,
      .phi = -0.175086,
      .Omega = 7.77823e15,
      .Gamma = 2.0376e15
    },
    .params[1] = {
      .bigA = 0.193194,
      .phi = -1.29864,
      .Omega = 6.17219e15,
      .Gamma = 0.409247e15
    },
  },
};

#endif
