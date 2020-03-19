#ifndef FDTD_FIELD_PROTOTYPES
#define FDTD_FIELD_PROTOTYPES

#include "fdtd_grid.h"

void EFieldUpdate(struct Grid *g);
void HFieldUpdate(struct Grid *g, int n);
void JFieldUpdate(struct Grid *g);
void DFTUpdate(struct Grid *g, int n);

void WriteDFTFile(struct Grid *g);
void NormalizeDFT(struct Grid *g);
void finishDFT(struct Grid *g);


#endif
