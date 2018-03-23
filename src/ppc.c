#include <R.h>
#include <Rdefines.h>

#include "ppc.h"
#include "ncdf.h"

/*  Main Program */
SEXP ppc(SEXP NCDF, SEXP Window, SEXP MassLimits, SEXP MinInt)
{
	ncdf_t *nc = new_ncdf(NCDF);
	matrix_t *mat = get_intensity_mat(nc);

	int min, max;
	if(!isNull(MassLimits)) {
		int *tmp = INTEGER(AS_INTEGER(MassLimits));
		min = tmp[0];
		max = tmp[1];
	} else {
		min = mat->mzmin;
		max = mat->mzmax;
	}

	int win     = INTEGER_VALUE(Window);
	int min_int = INTEGER_VALUE(MinInt);
	
	/* Rows = time, columns = m/z */
	int N = max - min + 1;
	SEXP MaxIntMatrix = PROTECT(allocMatrix(INTSXP, nc->nscans, N));
	int * maxm = INTEGER_POINTER(MaxIntMatrix);
	int * ans  = Calloc(nc->nscans, int);

	/* Looks for peaks for every mass (column) */
	for(int mz = min, i = 0; i < N; i++, mz++) {
		if(mz < mat->mzmin || mz > mat->mzmax)
			continue;
		int *x = mat->x + (mz - mat->mzmin) * mat->nr;
		peaks(x, win, nc->nscans, ans);
		/* assign maximum values to matrix */
		for(int j = 0; j < nc->nscans; j++) {
			maxm[j + i*nc->nscans] =
				(ans[j] == 1) && (x[j] >= min_int) ? x[j] : 0;
			ans[j] = 0;
		}
	}
	Free(ans);
	Free(nc);
	free_matrix(mat);
	UNPROTECT(1);
	return MaxIntMatrix;
}

/* This function implements PPC algorithm 
 * Arguments: x - Intensity vector
 *            ispan - Window
 *            n - length of x
 *            ans - output: 1 if it is a peak, 0 if not
 */            
int peaks(int *x, int ispan, int n, int *ans)
{
        int i, j;

        for(i = 0; i < ispan; i++)
                ans[i] = 0;

        for(i = n-ispan; i < n; i++)
                ans[i] = 0;

        i = ispan;

        while (i < n-ispan) {
                ans[i] = 1;
                j = i-ispan;
                while( ans[i] == 1 && j <= i+ispan) {
                        if(x[j] > x[i])
                                ans[i] = 0;
                        j++;
                        if(j == i)
                                j++;

                }
                if(ans[i] == 0)
                        i++;
                if(ans[i] == 1)
                        i += ispan-1;
        }

        return 1;
}
