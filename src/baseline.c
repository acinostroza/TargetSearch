/*
 * compute the baseline of a signal using quantiles
 */

#include "utils.h"
#include <Rdefines.h>

void
bsline(double *x, double *t, double qntl, double win, int step, int n,
	double *ret)
{
	/* compute the quantiles */
	double st = (double) step;
	double *q = Calloc(n, double);

	int qlen = qntl_win(x, t, qntl, win, step, n, q);

	/* linear interpolation of each point */
	for(int k = 1; k < qlen; k++) {
		for(int i = (k - 1) * step; i < k * step; i++) {
			double kk = (double) k, ii = (double) i;
			ret[i] = q[k-1] + (q[k] - q[k-1]) *
				(ii - (kk-1) * st) / st;
		}
	}

	/* fill the end with the last quantile value */
	for(int i = (qlen - 1) * step; i < n; i++)
		ret[i] = q[ qlen - 1 ];

	Free(q);
}

/* R interface for baseline correction */
SEXP baseline(SEXP X, SEXP T, SEXP Quantile, SEXP Window, SEXP Step)
{
	double *x   = REAL(AS_NUMERIC(X));
	double *t   = REAL(AS_NUMERIC(T));
	int win     = INTEGER_VALUE(Window);
	int step    = INTEGER_VALUE(Step);
	double qntl = NUMERIC_VALUE(Quantile);
	int n = GET_LENGTH(X);

	SEXP res = PROTECT(NEW_NUMERIC(n));
	double *bl = REAL(res);
	bsline(x, t, qntl, win, step, n, bl);
	UNPROTECT(1);
	return res;
}
