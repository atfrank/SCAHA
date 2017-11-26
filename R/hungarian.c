#include <R.h>
#include "hungarian.h"

void hungarian(int *nin, double *x, Sint *p)
{
  run(x, nin[0]);
  int i;
  for(i = 0; i < n; i++)
	 p[i] = xy[i];
}
