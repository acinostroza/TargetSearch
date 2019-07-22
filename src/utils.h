#ifndef UTILS_H
#define UTILS_H

#include <R.h>

#ifndef USING_R
	#include <stdlib.h>
	#include <string.h>
	#include <math.h>
	#define Calloc(n,t) (t *) calloc( n, sizeof(t) )
	#define Free(x) free(x)
	#define Memzero(p,n) memset(p, 0, (n) * sizeof(*p))
	#pragma message "not using R"
#endif

void moving(const int *, int, int, double *);
void convolve(const int *, int, const double *, int, double *);
int find_peak_diff(const double *, int, int *);
void refine_peak(const int *, int, int, int *, int);
double * gaussian_coef(int);
int peak_detection_ppc(int *, int, int, int *);
int qntl_win(double *, double *, double, double, int, int, double *);

#endif
