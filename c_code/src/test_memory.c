#include <stdlib.h>
#include <stdio.h>
#include "fdtd_proto.h"
#include "fdtd_macro.h"

int main()
{

  int m,n;
  m = 500;
  n = 501;

  struct Grid *g;
  g = AllocateGridMemory();

  InitializeFdtd(g, 0);

/*  double **ptr1,**ptr2,**ptr3, **ptr4, **ptr5, **ptr6;
  ptr1 = AllocateMemory(m, n, 0.0);
  ptr2 = AllocateMemory(m, n, 0.0);
  ptr3 = AllocateMemory(m, n, 0.0);
  ptr4 = AllocateMemory(m, n, 0.0);
  ptr5 = AllocateMemory(m, n, 0.0);
  ptr6 = AllocateMemory(m, n, 0.0);

  ex = AllocateMemory(m, n, 1.0);
  ey = AllocateMemory(m, n, 1.0);
  hz = AllocateMemory(m, n, 1.0);
*/
  printf("ex[0][0]: %f\n", ex[0][0]);
  freeGrid(g);
  
  return(0);
}
