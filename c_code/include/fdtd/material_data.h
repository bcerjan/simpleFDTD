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
    then fit with a simple Drude model (simultaneously for Re and Im components):
      er[f] = epsInfinity - fp^2 / ( f^2 + i*g*f ) .

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
static const double materialData[5][6] = {
  //0, Al, Cheng 2016
  {4.01771, 3.77e+7, 1.0, 2.65e-8, 22.4842e+14,5.309e+14},
  //1, Gold, Johnson and Christy 1974
  {2.32842, 4.11e+7, 1.0, 2.44e-8, 12.0356e+14,2.3138e+14},
  //2, Silver, Wu XXX
  {7.19906e-9, 6.30e+7, 1.0, 1.59e-8, 9.71771e+14, 0.314692e+14},
  //3, Copper, Johnson and Christy 1974
  {1.81405e-7, 5.69e+7, 1.0, 1.68e-8, 17.9867e+14, 12.9089e+14},
  //4, Silica, fixed eps = 2.136 -> n = 1.46
  {2.136, 0.0, 1.0, 0.0, HUGE_VAL, 0.0}

};


#endif
