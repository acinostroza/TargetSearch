/* functions for baseline correction */

#include <R.h>
#include <Rinternals.h>

#include "hpf.h"

/* a first order high pass filter */
void hpf(double *x, double *y, int *n, double *a) {
	int i = 0;
	y[0] = x[0];
	for(i = 1; i < *n; i++)
		y[i] = *a * (y[i-1] + x[i] - x[i-1]);
}

/* windowing step*/

void windowing (int *x, int *idx, int *w, int *n, int *nidx) {
        int i,j,a,b;
        for(i = 0; i < *nidx; i++) {
                a = ( idx[i] - *w < 1 )  ?  1 : idx[i] - *w;
                b = ( idx[i] + *w > *n ) ? *n : idx[i] + *w;
                for(j = a; j <= b; j++)
                        x[j - 1] = 1;
        }
}
