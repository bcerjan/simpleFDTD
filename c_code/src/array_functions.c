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

#include <math.h>

/* File containing array functions (max, min) */

double ArrayMax(double **ptr, int xWidth, int yWidth)
{
  int i,j;
  double max = ptr[0][0];
  for (i = 0; i < xWidth; i++) {
    for (j = 0; j < yWidth; j++) {
      if (ptr[i][j] > max) {
        max = ptr[i][j];
      } /* ifBlock */
    } /* jForLoop */
  } /* iForLoop */
  return max;
}

double ArrayMin(double **ptr, int xWidth, int yWidth)
{
  int i,j;
  double min = ptr[0][0];
  for (i = 0; i < xWidth; i++) {
    for (j = 0; j < yWidth; j++) {
      if (ptr[i][j] < min) {
        min = ptr[i][j];
      } /* ifBlock */
    } /* jForLoop */
  } /* iForLoop */
  return min;
}

double AbsArrayMax(double **ptr, int xWidth, int yWidth)
{
  int i,j;
  double max = ptr[0][0];
  double temp = 0.0;
  for (i = 0; i < xWidth; i++) {
    for (j = 0; j < yWidth; j++) {
      temp = fabs(ptr[i][j]);
      if (temp > max) {
        max = temp;
      } /* ifBlock */
    } /* jForLoop */
  } /* iForLoop */
  return max;
}

double AbsArrayMin(double **ptr, int xWidth, int yWidth)
{
  int i,j;
  double min = ptr[0][0];
  double temp = 0.0;
  for (i = 0; i < xWidth; i++) {
    for (j = 0; j < yWidth; j++) {
      temp = fabs(ptr[i][j]);
      if (temp < min) {
        min = temp;
      } /* ifBlock */
    } /* jForLoop */
  } /* iForLoop */
  return min;
}
