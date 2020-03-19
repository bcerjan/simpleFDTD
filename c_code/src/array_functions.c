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
