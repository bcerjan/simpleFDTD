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

#ifndef OUTPUTFUNCS
#define OUTPUTFUNCS

struct AuxIndexFields {
  // Storage for fields(n,w):
  double **nTranDFT;
  double **nReflDFT;
  double **nReEy;
  double **nImEy;
  double **nReHz;
  double **nImHz;

  // Number of indices we're using:
  int numIndex;
  // Number time steps (in case we need it later for verification)
  int maxSteps;

  int dftFreqs;
};

void WriteDFTFile(struct AuxIndexFields *Fields);

#endif
