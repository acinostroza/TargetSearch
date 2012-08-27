/* Function to correct netCDF data when mass values are not integers as
 * expected by TargetSearch. 'cdffix_core' assign each mass value to the
 * nearest integer taking the highest intensity in case of ambiguities.
 * There is a limit in the number of ambiguities give by constant MAX_ASSIGNED
 * (currently 3).

 * Note:
     This function is not intended for high mass accuracy data.
*/

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#define  MAX_ASSIGNED  3
#define  TRUE          1
#define  FALSE         0

/* fixed cdf data*/
struct cdf {
	int *mass;
	int *in;
	int *scan_idx;
	int *p_count;
	int np;       /* no. points */
	int ns;       /* no. scans  */
};

struct cdf *
init_cdf(int ns, int np)
{
	struct cdf *x = Calloc(1, struct cdf);

	x->mass     = Calloc(np, int);
	x->in       = Calloc(np, int);
	x->scan_idx = Calloc(ns, int);
	x->p_count  = Calloc(ns, int);
	x->ns       = ns;
	x->np       = np;
	return x;
}

void
free_cdf(struct cdf *x)
{
	Free(x->mass);
	Free(x->in);
	Free(x->scan_idx);
	Free(x->p_count);
	Free(x);
}

int
cdffix_core(struct cdf *x, double *m, int *in, int *scan_idx, int *p_count)
{
	int i, j, k, p=0, count = 0;
	int tmp;

	for(i = 0; i < x->ns; i++) {
		count = x->p_count[i] = 0;
		for(j = 0; j < p_count[i]; j++) {
			k = j + scan_idx[i];
			tmp = (int) fround(m[k], 0.0);

			if(p == 0) {
				x->mass[p] = tmp;
				x->in[p]   = in[k];
				x->p_count[i]++;
				p++;
				continue;
			}

			if(tmp == x->mass[p-1]) {
				count++;
				/* abort function when there're too many
				 * masses assigned to the same integer */
				if(count >= MAX_ASSIGNED)
					return FALSE;

				if(in[k] > x->in[p-1])
					x->in[p-1] = in[k];
			} else {
				x->mass[p] = tmp;
				x->in[p]   = in[k];
				x->p_count[i]++;
				count = 0;
				p++;
			}
		}
	}
	x->np = p;
	x->scan_idx[0] = 0;
	for(i = 1; i < x->ns; i++)
		x->scan_idx[i] = x->scan_idx[i-1] + x->p_count[i-1];
	return TRUE;
}

SEXP
cdffix(SEXP mass, SEXP intensity, SEXP scanIndex, SEXP pointCount)
{
	double *m;
	int    *in, *scan_idx, *p_count;
	int    ns, np;
	struct cdf *x;
	SEXP   M, I, S, P;
	SEXP   res;

	mass       = AS_NUMERIC(mass);
	intensity  = AS_INTEGER(intensity);
	scanIndex  = AS_INTEGER(scanIndex);
	pointCount = AS_INTEGER(pointCount);

	m        = NUMERIC_POINTER(mass);
	in       = INTEGER_POINTER(intensity);
	scan_idx = INTEGER_POINTER(scanIndex);
	p_count  = INTEGER_POINTER(pointCount);

	ns = GET_LENGTH(scanIndex);
	np = GET_LENGTH(mass);

	x = init_cdf(ns, np);
	if(cdffix_core(x, m, in, scan_idx, p_count) != TRUE)
		return R_NilValue;

	PROTECT(res = allocVector(VECSXP, 4));
	PROTECT(M   = NEW_INTEGER(x->np));
	PROTECT(I   = NEW_INTEGER(x->np));
	PROTECT(S   = NEW_INTEGER(x->ns));
	PROTECT(P   = NEW_INTEGER(x->ns));

	Memcpy(INTEGER_POINTER(M), x->mass,     x->np);
	Memcpy(INTEGER_POINTER(I), x->in,       x->np);
	Memcpy(INTEGER_POINTER(S), x->scan_idx, x->ns);
	Memcpy(INTEGER_POINTER(P), x->p_count,  x->ns);

	SET_VECTOR_ELT(res, 0, M);
	SET_VECTOR_ELT(res, 1, I);
	SET_VECTOR_ELT(res, 2, S);
	SET_VECTOR_ELT(res, 3, P);

	UNPROTECT(5);
	free_cdf(x);
	return res;
}

