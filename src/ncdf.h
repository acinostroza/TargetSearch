#ifndef _NCDF_H
#define _NCDF_H

#include <R.h>
#include <Rdefines.h>

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

typedef struct
{
	int *x;
	int nc;
	int nr;
	int mzmin;
	int mzmax;
	int alloc;    /* total memory allocated or 0 */
} matrix_t;

matrix_t * get_intensity_mat(ncdf_t *);
ncdf_t   * new_ncdf(SEXP);
matrix_t * from_matrix(SEXP);
void free_matrix(matrix_t *);

/* fix a CDF with non-integer mass values (aka non-nominal mass) */
SEXP cdffix(SEXP, SEXP);

/* get the intensity matrix of a NCDF object */
SEXP ncdfToMatrix(SEXP NCDF, SEXP massRange);

#endif
