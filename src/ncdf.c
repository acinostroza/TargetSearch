#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "ncdf.h"

#define  TRUE          1
#define  FALSE         0

/* get an element from a list from their names
 * adapted from "writing R extensions"
 */
static SEXP get(SEXP list, const char *str)
{
	SEXP elem = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	for(int i = 0; i < GET_LENGTH(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elem = VECTOR_ELT(list, i);
			break;
		}
	return elem;
}

/* functions to manipulate CDF data in form of a list into a C structure */

/* returns the mz limits */
static void mass_range(ncdf_t *x, int *min, int *max)
{
	int *mz = x->mass;
	*min = *max = mz[0];
	for(int i = 1; i < x->npoints; i++)
	{
		if(mz[i] > *max)
			*max = mz[i];
		if(mz[i] < *min)
			*min = mz[i];
	}
}

/* matrix of intensities */
matrix_t * get_intensity_mat(ncdf_t *x)
{
	int min, max, n, *z;
	matrix_t *mat = R_Calloc(1, matrix_t);
	mass_range(x, &min, &max);
	mat->mzmin = min;
	mat->mzmax = max;
	mat->nc    = (max - min + 1);
	mat->nr    = x->nscans;
	mat->alloc = n = (max - min + 1) * x->nscans;
	z = R_Calloc(n, int);

	for(int s = 0; s < x->nscans; s++) {
		for(int i = 0; i < x->p_count[s]; i++) {
			z[(x->mass[x->scan_idx[s] + i]-min) * x->nscans + s] = x->in[x->scan_idx[s] + i];
		}
	}
	mat->x = z;
	return mat;
}

/* transform a R matrix of intensities into a matrix_t struct */
matrix_t * from_matrix(SEXP Matrix)
{
	SEXP dim = GET_DIM(Matrix);
	matrix_t *ret = NULL;

	if(isNull(dim))
		return ret;
	ret = new_mat_alloc(INTEGER(dim)[1], INTEGER(dim)[0], INTEGER(AS_INTEGER(Matrix)));
	return ret;
}

/* transform R object into a struct ncdf_t. It expects that the m/z and intensity
 * values have been already rounded in R (as it is not necessary to implement this
 * in C).
 * To free the memory, call free_ncdf
 */
ncdf_t * new_ncdf(SEXP NCDF)
{
	ncdf_t *cdf = R_Calloc(1, ncdf_t);
	SEXP RT         = get(NCDF, "rt");
	SEXP RI         = get(NCDF, "ri");
	SEXP PointCount = get(NCDF, "point_count");
	SEXP ScanIndex  = get(NCDF, "scanindex");
	SEXP MZ         = get(NCDF, "mz");
	SEXP Intensity  = get(NCDF, "intensity");
	cdf->nscans     = GET_LENGTH(RT);
	cdf->npoints    = GET_LENGTH(MZ);

	/* retention index might be null if not set */
	cdf->ri = isNull(RI) ? NULL : REAL(AS_NUMERIC(RI));
	cdf->rt = REAL(AS_NUMERIC(RT));
	cdf->p_count = INTEGER(AS_INTEGER(PointCount));
	cdf->scan_idx = INTEGER(AS_INTEGER(ScanIndex));
	cdf->mass = INTEGER(AS_INTEGER(MZ));
	cdf->in = INTEGER(AS_INTEGER(Intensity));
	cdf->alloc = 0; /* flag to indicate allocated memory */
	return cdf;
}

/******************************************************************************
 * Functions to transform netCDF data to nominal mass when the mass values are
 * not integers. 'cdf_to_nominal' assign each mass value to the nearest integer
 * summing up the intensities of duplicated mass values.
 * Note:
 * These functions are intended for centroid data, though it may also work
 * for profile data. There is of course loss of information due to the rounding.
 ******************************************************************************/

/* allocate a new ncdf_t object using transient allocation methods. */
static ncdf_t *
alloc_cdf(int ns, int np)
{
	ncdf_t *x   = R_Calloc(1, ncdf_t);

	x->mass     = R_Calloc(np, int);
	x->in       = R_Calloc(np, int);
	x->scan_idx = R_Calloc(ns, int);
	x->p_count  = R_Calloc(ns, int);
	x->rt       = R_Calloc(ns, double);
	x->ri       = R_Calloc(ns, double);
	x->nscans   = ns;
	x->npoints  = np;
	x->alloc    = 1;
	return x;
}

/* free allocated ncdf_t object. checks whether transiently allocation memory
 * was used, ie, does not free shared memory.
 */
void
free_ncdf(ncdf_t *x)
{
	if(x == NULL)
		return;
	if(x->alloc == 1) {
		R_Free(x->rt);
		R_Free(x->ri);
		R_Free(x->mass);
		R_Free(x->in);
		R_Free(x->scan_idx);
		R_Free(x->p_count);
	}
	R_Free(x);
}

/* main function to transform CDF data to nominal mass. Duplicated intensities are
 * summed up */
static int
nominal_core(ncdf_t *dest, ncdf_t *src)
{
	int k = 0;

	for(int i = 0; i < src->nscans; i++) {
		int si = src->scan_idx[i];
		int pc = src->p_count[i];

		for(int j = si; j < si + pc; j++) {
			if(j == si || src->mass[j] != src->mass[j - 1]) {
				dest->mass[k] = src->mass[j];
				dest->in[k] = src->in[j];
				dest->p_count[i]++;
				k++;
			}
			else {
				dest->in[k - 1] += src->in[j];
			}
		}
	}
	dest->npoints = k;
	dest->scan_idx[0] = 0;

	for(int i = 0; i < dest->nscans; i++) {
		if(i > 0)
			dest->scan_idx[i] = dest->scan_idx[i-1] + dest->p_count[i-1];
		dest->ri[i] = (src->ri != NULL) ? src->ri[i] : 0.0;
		dest->rt[i] = src->rt[i];
	}
	return TRUE;
}

/* takes a ncdf_t object and returns a NCDF list */
SEXP ncdf_sexp(ncdf_t *x)
{
	SEXP res   = PROTECT(allocVector(VECSXP, 6));
	SEXP names = PROTECT(allocVector(STRSXP, 6));

	SET_VECTOR_ELT(res, 0, NEW_NUMERIC(x->nscans));
	SET_VECTOR_ELT(res, 1, NEW_NUMERIC(x->nscans));
	SET_VECTOR_ELT(res, 2, NEW_INTEGER(x->nscans));
	SET_VECTOR_ELT(res, 3, NEW_INTEGER(x->nscans));
	SET_VECTOR_ELT(res, 4, NEW_INTEGER(x->npoints));
	SET_VECTOR_ELT(res, 5, NEW_INTEGER(x->npoints));

	Memcpy(REAL(VECTOR_ELT(res, 0)), x->rt, x->nscans);
	Memcpy(REAL(VECTOR_ELT(res, 1)), x->ri, x->nscans);
	Memcpy(INTEGER(VECTOR_ELT(res, 2)), x->scan_idx, x->nscans);
	Memcpy(INTEGER(VECTOR_ELT(res, 3)), x->p_count, x->nscans);
	Memcpy(INTEGER(VECTOR_ELT(res, 4)), x->mass, x->npoints);
	Memcpy(INTEGER(VECTOR_ELT(res, 5)), x->in,   x->npoints);

	SET_STRING_ELT(names, 0, mkChar("rt"));
	SET_STRING_ELT(names, 1, mkChar("ri"));
	SET_STRING_ELT(names, 2, mkChar("scanindex"));
	SET_STRING_ELT(names, 3, mkChar("point_count"));
	SET_STRING_ELT(names, 4, mkChar("mz"));
	SET_STRING_ELT(names, 5, mkChar("intensity"));
	setAttrib(res, R_NamesSymbol, names);
	return res;
}

/*************************/
/* .Call interface to R  */
/*************************/

/* Transform a CDF to nominal mass, ie, non-integer mass values
 * Args:
 *   NCDF: a list of named elements holding the CDF structure:
 *      rt : Retention time (double)
 *      ri : Retention time index (double)
 *      mz: m/z values (integer, duplicated values allowed)
 *      intensity: intensity values (integer)
 *      scanindex: scan index (integer)
 *      point_count: number of points per scan (integer).
 *   MA: max assigned to the same integer (to detect high mass accuracy)
 * Output:
 *   A similar list as the input, but with integer m/z values
 * Note:
 *   The m/z value must be rounded before correcting the cdf data
 *   (It would look like the m/z values are duplicated because of the
 *   rounding).
 */

SEXP
nominal(SEXP NCDF)
{
	ncdf_t *src = new_ncdf(NCDF);
	ncdf_t *dst = alloc_cdf(src->nscans, src->npoints);
	SEXP res = R_NilValue;

	if(nominal_core(dst, src) == TRUE)
		res = ncdf_sexp(dst);

	free_ncdf(dst);
	free_ncdf(src);

	if(!isNull(res))
		UNPROTECT(2); /* unprotects call to ncdf_sexp */

	return res;
}

/* takes a CDF object and return a matrix *
 * this function replaces peakCDFextraction */
SEXP ncdf_to_matrix(SEXP NCDF, SEXP massRange)
{
	ncdf_t *nc = new_ncdf(NCDF);
	matrix_t *mat = get_intensity_mat(nc);
	int *tmp = INTEGER(AS_INTEGER(massRange));
	int min = tmp[0], max = tmp[1];
	SEXP ansMat = PROTECT(allocMatrix(INTSXP, nc->nscans, max - min + 1));
	int *z = INTEGER_POINTER(ansMat);
	for(int mz = min, i = 0; mz <= max; i++, mz++) {
		if(mz < mat->mzmin || mz > mat->mzmax)
			continue;
		int *x = mat->x + (mz - mat->mzmin) * mat->nr;
		for(int j = 0; j < nc->nscans; j++)
			z[j + i*nc->nscans] = x[j];
	}
	R_Free(nc);
	free_mat(mat);
	UNPROTECT(1);
	return ansMat;
}
