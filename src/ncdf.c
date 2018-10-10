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
	matrix_t *mat = Calloc(1, matrix_t);
	mass_range(x, &min, &max);
	mat->mzmin = min;
	mat->mzmax = max;
	mat->nc    = (max - min + 1);
	mat->nr    = x->nscans;
	mat->alloc = n = (max - min + 1) * x->nscans;
	z = Calloc(n, int);

	for(int s = 0; s < x->nscans; s++) {
		for(int i = 0; i < x->p_count[s]; i++) {
			z[(x->mass[x->scan_idx[s] + i]-min) * x->nscans + s] = x->in[x->scan_idx[s] + i];
		}
	}
	mat->x = z;
	return mat;
}

/* tranform a R matrix of intensities into a matrix_t struct */
matrix_t * from_matrix(SEXP Matrix)
{
	SEXP dim = GET_DIM(Matrix);
	matrix_t *ret = NULL;

	if(isNull(dim))
		return ret;
	ret = new_mat_alloc(INTEGER(dim)[1], INTEGER(dim)[0], INTEGER(AS_INTEGER(Matrix)));
	return ret;
}

/* transform R object into a struct ncdf_t. It expects an already fixed cdf list
 * in which masses and intensities are integers. remember to call `Free' on the
 * generated object
 */
ncdf_t * new_ncdf(SEXP NCDF)
{
	ncdf_t *cdf = Calloc(1, ncdf_t);
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
 * Functions to correct netCDF data when mass values are not integers as
 * expected by TargetSearch. 'cdffix_core' assign each mass value to the
 * nearest integer taking the highest intensity in case of ambiguities.
 * There is a limit in the number of ambiguities give by constant MAX_ASSIGNED
 * (currently 3).
 *
 * Note:
 * This functions are not intended for high mass accuracy data.
 ******************************************************************************/

/* allocate a new ncdf_t object using transient allocation methods.
 * Note that 'new_ncdf' uses shared memory */
static ncdf_t *
alloc_cdf(int ns, int np)
{
	ncdf_t *x   = Calloc(1, ncdf_t);

	x->mass     = Calloc(np, int);
	x->in       = Calloc(np, int);
	x->scan_idx = Calloc(ns, int);
	x->p_count  = Calloc(ns, int);
	x->rt       = Calloc(ns, double);
	x->ri       = Calloc(ns, double);
	x->nscans   = ns;
	x->npoints  = np;
	x->alloc    = 1;
	return x;
}

/* free allocated ncdf_t object. To be used in an transiently allocation object */
void
free_cdf(ncdf_t *x)
{
	if(x == NULL)
		return;
	if(x->alloc == 1) {
		Free(x->rt);
		Free(x->ri);
		Free(x->mass);
		Free(x->in);
		Free(x->scan_idx);
		Free(x->p_count);
	}
	Free(x);
}

/* main function to fix a cdf data */
static int
cdffix_core(ncdf_t *dest, ncdf_t *src, int max_assigned)
{
	int i, j, k, p=0, count = 0;
	int *scan_idx = src->scan_idx, *p_count = src->p_count;
	int *in = src->in, *m = src->mass;
	ncdf_t *x = dest;

	for(i = 0; i < x->nscans; i++) {
		count = x->p_count[i] = 0;
		for(j = 0; j < p_count[i]; j++) {
			k = j + scan_idx[i];

			if(p == 0) {
				x->mass[p] = m[k];
				x->in[p]   = in[k];
				x->p_count[i]++;
				p++;
				continue;
			}

			if(m[k] == x->mass[p-1]) {
				count++;
				/* abort function when there're too many
				 * masses assigned to the same integer */
				if(count >= max_assigned)
					return FALSE;

				if(in[k] > x->in[p-1])
					x->in[p-1] = in[k];
			} else {
				x->mass[p] = m[k];
				x->in[p]   = in[k];
				x->p_count[i]++;
				count = 0;
				p++;
			}
		}
	}
	x->npoints = p;

	x->scan_idx[0] = 0;
	for(i = 0; i < x->nscans; i++) {
		if(i > 0)
			x->scan_idx[i] = x->scan_idx[i-1] + x->p_count[i-1];
		x->ri[i] = (src->ri != NULL) ? src->ri[i] : 0.0;
		x->rt[i] = src->rt[i];
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

/* Fix a CDF file with non-integer mass values (no high mass accuracy)
 * Args:
 *   NCDF: a list of named elements holding the CDF structure:
 *      rt : Retention time (double)
 *      ri : Retention time index (double)
 *      mz: m/z values (integer)
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
cdffix(SEXP NCDF, SEXP MA)
{
	ncdf_t *dst, *nc = new_ncdf(NCDF);
	SEXP  res  = R_NilValue;
	dst = alloc_cdf(nc->nscans, nc->npoints);
	int max_assigned = isNull(MA) ? 10000000 : INTEGER(MA)[0];

	if(cdffix_core(dst, nc, max_assigned) == TRUE)
		res = ncdf_sexp(dst);
	free_cdf(dst);

	Free(nc);
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
	Free(nc);
	free_mat(mat);
	UNPROTECT(1);
	return ansMat;
}
