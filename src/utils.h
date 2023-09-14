#ifndef UTILS_H
#define UTILS_H

#include <R.h>

void moving(const int *, int, int, double *);
void convolve(const int *, int, const double *, int, double *);
int find_peak_diff(const double *, int, int *);
void refine_peak(const int *, int, int, int *, int);
double * gaussian_coef(int);
int peak_detection_ppc(int *, int, int, int *);
int qntl_win(double *, double *, double, double, int, int, double *);
int binsearch(double *, double, int);

#endif
