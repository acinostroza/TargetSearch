#ifndef _NCDF_H
#define _NCDF_H

#include <R.h>
#include <Rdefines.h>
#include "matrix.h"

typedef struct
{
	double *rt;     /* Retention Time */
	double *ri;     /* Retention Index */
	int *p_count;   /* Point Count */
	int *scan_idx;  /* Scan Index */
	int *mass;      /* m/z */
	int *in;        /* intensity */
	int nscans;
	int npoints;
} ncdf_t;

matrix_t * get_intensity_mat(ncdf_t *);
ncdf_t   * new_ncdf(SEXP);
matrix_t * from_matrix(SEXP);

/* fix a CDF with non-integer mass values (aka non-nominal mass) */
SEXP cdffix(SEXP, SEXP);

/* get the intensity matrix of a NCDF object */
SEXP ncdf_to_matrix(SEXP NCDF, SEXP massRange);

#endif
