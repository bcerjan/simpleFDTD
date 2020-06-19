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

#include <stdio.h>
#include "output_funcs.h"

// Function to put the licensing info at the top of auto-generated header files
// Assumes a stream has already been opened (you should have had to do to
// generate the pointer anyway)
void writeLicenseHeader(FILE *filePtr) {
  fprintf(filePtr, "/**\n \
    Copyright (c) 2020 Ben Cerjan\n\n \
    This file is part of simpleFDTD.\n\n \
    simpleFDTD is free software: you can redistribute it and/or modify\n \
    it under the terms of the GNU Affero General Public License as published by\n \
    the Free Software Foundation, either version 3 of the License, or\n \
    (at your option) any later version.\n\n \
    simpleFDTD is distributed in the hope that it will be useful,\n \
    but WITHOUT ANY WARRANTY; without even the implied warranty of\n \
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n \
    GNU Affero General Public License for more details.\n\n \
    You should have received a copy of the GNU Affero General Public License\n \
    along with simpleFDTD.  If not, see <https://www.gnu.org/licenses/>.\n \
    **/\n\n");

    return;

}


// Function to write 2-d field to header -- does not do initial setup, just
// outputs rows for an initialized array
void writeFieldVals (FILE *filePtr, double **fieldPtr, int iSteps, int jSteps) {
  int i,j;

  for (i = 0; i < iSteps; i++){
    for (j = 0; j < jSteps; j++) {
      fprintf(filePtr, "%.17g,\n", fieldPtr[i][j]);
    } /* jForLoop */
    if ( i < iSteps - 1 ) {
      fprintf(filePtr, "},\n{");
    } else {
      fprintf(filePtr, "}};\n" );
    } /* if */
  } /* iForLoop */

  return;
}

/*
 Function to write header files for our DFT normalization
 For the transmitted wave, we only need to store P(index, w), but for the
 reflected wave, we need to store the complex Ey(index,w) and Hz(index,w)
 for later analysis
*/
void WriteDFTFile (struct AuxIndexFields *Fields) {
  int i,j;
  char reflFilename[100] = "../include/fdtd/empty_refl_data.h";
  char tranFilename[100] = "../include/fdtd/empty_tran_data.h";
  FILE *reflDataPtr, *tranDataPtr;
  int numInd = Fields->numIndex;
  int numFreqs = Fields->dftFreqs;

  // Write to header file for use later
  reflDataPtr = fopen(reflFilename, "w");
  tranDataPtr = fopen(tranFilename, "w");

  // Step 0: Add licensing text:
  writeLicenseHeader(reflDataPtr);
  writeLicenseHeader(tranDataPtr);

  // First add blocking definitions (just in case...)
  fprintf(reflDataPtr, "#ifndef REFL_EMPTY_DATA\n#define REFL_EMPTY_DATA\n");
  fprintf(tranDataPtr, "#ifndef TRAN_EMPTY_DATA\n#define TRAN_EMPTY_DATA\n");

  // Second add lines for number of freqs and time steps:
  fprintf(reflDataPtr, "#define reflSteps %i\n", Fields->maxSteps);
  fprintf(reflDataPtr, "#define reflFreqs %i\n", numFreqs);
  fprintf(tranDataPtr, "#define tranSteps %i\n", Fields->maxSteps);
  fprintf(tranDataPtr, "#define tranFreqs %i\n", numFreqs);

  // Do the Transmittance and Reflected Flux arrays first because they're simpler:
  fprintf(tranDataPtr, "static const double emptyTranDFT[%i][%i] = {\n{\n", numInd, numFreqs );
  fprintf(reflDataPtr, "static const double emptyReflDFT[%i][%i] = {\n{\n", numInd, numFreqs );

  // Write in the values:
  writeFieldVals(reflDataPtr, Fields->nReflDFT, numInd, numFreqs);
  writeFieldVals(tranDataPtr, Fields->nTranDFT, numInd, numFreqs);

  // And close the file, since we're now done with it
  fprintf(tranDataPtr, "\n" );
  fprintf(tranDataPtr, "#endif"); // End if block
  fclose(tranDataPtr);

  /* For Reflectance data, we need four more arrays, two for Re/Im(Ey(index,w,{x})) and
     two for Re/Im(Hz(index,w,{x})). The position dependence is suppressed, as they
     are actually constant with position. Note that this assumption assumes that
     the source is a uniform line source. If it is not, this does not work and
     needs to retain the position dependence.
   */
   
  fprintf(reflDataPtr, "static const double emptyReEyRefl[%i][%i] = {\n{\n", numInd, numFreqs );
  writeFieldVals(reflDataPtr, Fields->nReEy, numInd, numFreqs);

  // End this array, start next one, Im(Ey):
  fprintf(reflDataPtr, "\nstatic const double emptyImEyRefl[%i][%i] = {\n{\n", numInd, numFreqs );
  writeFieldVals(reflDataPtr, Fields->nImEy, numInd, numFreqs);

  // End this array, start next one, Re(Hz):
  fprintf(reflDataPtr, "\nstatic const double emptyReHzRefl[%i][%i] = {\n{\n", numInd, numFreqs );
  writeFieldVals(reflDataPtr, Fields->nReHz, numInd, numFreqs);

  // End this array, start next one, Im(Hz):
  fprintf(reflDataPtr, "\nstatic const double emptyImHzRefl[%i][%i] = {\n{\n", numInd, numFreqs );
  writeFieldVals(reflDataPtr, Fields->nImHz, numInd, numFreqs);

  // End the array:
  fprintf(reflDataPtr, "\n" );
  fprintf(reflDataPtr, "#endif"); // End if block
  fclose(reflDataPtr);

  return;
}
