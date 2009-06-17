#include <R.h>
#include <Rdefines.h>

int peaks(int *, int, int, int *);

/*  Main Program */

SEXP ppc(SEXP MassValues, SEXP IntensityValues, SEXP PointCount,
	SEXP ScanIndex, SEXP Window, SEXP MassLimits, SEXP MinInt) {
	
	int min, max, win, N, min_int;
	int scan_number, *scan_index, *point_count, *intensity, *mass;
	
	int *ans, *maxm; /* pointer to max intensity matrix */

	int s, i;
	
	SEXP MaxIntMatrix;

	/* R objects created in the C code have to be reported using the PROTECT 
	 * macro on a pointer to the object. This tells R that the object is in 
	 * use so it is not destroyed. */
	PROTECT(MassValues = AS_INTEGER(MassValues));
	PROTECT(IntensityValues = AS_INTEGER(IntensityValues));
	PROTECT(PointCount = AS_INTEGER(PointCount));
	PROTECT(ScanIndex = AS_INTEGER(ScanIndex));
	PROTECT(MassLimits = AS_INTEGER(MassLimits));
		
	min = INTEGER_POINTER(MassLimits)[0];
	max = INTEGER_POINTER(MassLimits)[1];
	win = INTEGER_VALUE(Window);
	min_int = INTEGER_VALUE(MinInt);
	
	N = max - min + 1;
	
	scan_number = GET_LENGTH(ScanIndex);
	
	scan_index  = INTEGER_POINTER(ScanIndex);
	point_count = INTEGER_POINTER(PointCount);
	intensity   = INTEGER_POINTER(IntensityValues);
	mass        = INTEGER_POINTER(MassValues);

	ans  = Calloc(N * scan_number, int);

	/* Rows = time, columns = m/z */
	PROTECT(MaxIntMatrix = allocMatrix(INTSXP, scan_number, N));
	maxm = INTEGER_POINTER(MaxIntMatrix);

	/* fill Intensity matrix with raw data */
	for(s = 0; s < scan_number; s++) {
		/* Make sure that the matrix is filled with zeros */
		for (i = 0; i < N; i++)
			maxm[i*scan_number+s] = 0;

		for (i = 0; i < point_count[s]; i++)
			 if(mass[scan_index[s]+i] >= min && mass[scan_index[s]+i] <= max)
				maxm[(mass[scan_index[s] + i]-min) * scan_number + s] = intensity[scan_index[s] + i];
	}
	/* Looks for peaks for every mass (column) */	
	for(i = 0; i < N; i++)
		peaks(maxm + i*scan_number, win, scan_number, ans + i*scan_number);

	/* Set everything that is not a peak to zero.
	* Do not look for points that we already know are zero, ie, absent in
	* raw data */
	for(s = 0; s < scan_number; s++)
		for (i = 0; i < point_count[s]; i++)
			if((ans[(mass[scan_index[s] + i]-min) * scan_number + s] == 0) ||
				(intensity[scan_index[s] + i] < min_int))
				maxm[(mass[scan_index[s] + i]-min) * scan_number + s] = 0;

	UNPROTECT(6);
	Free(ans);
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
