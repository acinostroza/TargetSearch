/* Change Log:
	- 2008.08.18: new algorithm to calculate the differences of smooth intensities. This fixed a integer
		overflow bug when a big smoothing window is used.
	- 2008.08.13: A new condition was included to detect peaks: if the difference at time T is zero,
		the differences at times T-1 and T+1 are also checked for change of sign (from positive to negative)
	- 2008.06.09: Fixed. MassLimits were not properly checked.
*/

#include <R.h>
#include <Rdefines.h>

/* local maximun window size (plus minus) */
#define LM_WIN 4

/* Function prototypes */

int get_int(int , int *, int, int, int *, int *,int *, int *);
void sum(int *x, int *t, int n);
int check(int *, int, int);
int get_max_int(int, int, int, int *, int *, int *, int *, int *);

/*  Main Program */

SEXP peak_finding(SEXP MassValues, SEXP IntensityValues, SEXP PointCount,
	SEXP ScanIndex, SEXP Window, SEXP MassLimits, SEXP MinInt) {
	
	int min, max, win, N, min_int;
	int scan_number, *scan_index, *point_count, *intensity, *mass;
	
	int *diff, *maxm; /* pointer to max intensity matrix */

	int s, i, max_int, s_max_int;
	int *x, *y;
	
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

	diff = (int *) R_alloc( N * (scan_number - 1), sizeof(int));
	x    = (int *) R_alloc( N, sizeof(int));
	y    = (int *) R_alloc( N, sizeof(int));

	PROTECT(MaxIntMatrix = allocMatrix(INTSXP, scan_number, N));
	maxm = INTEGER_POINTER(MaxIntMatrix);

	/* Initialize intensity matrix */
	for(s = 0; s < scan_number*N; s++)
		maxm[s] = 0;
	
	/* The differences of intensities for each mass are computed.
	 * Xmean(s,w) = { sum_(k=s-w)^(s+w) X(k) } / { 2w+1 } , where X(s) is the intensity at scan s
	 * D(s,w) = Xmean(s,w) - Xmean(s-1,w) = {X(s+w) - X(s-w-1)} / { 2w+1 }
	* Note that is not necessary to calculate the average of the intensities
	*/
	
	for(s = 0; s < scan_number - 1; s++) {

		for (i = 0; i < N; i++)
			x[i] = y[i] = 0;
	
		/* Get scan s+win  */
		if ( s+win < scan_number)
			get_int(s+win, y, min, max, mass, intensity, scan_index, point_count);
		
		/* Get scan s-win-1 */
		if (s-win-1 >= 0)
			get_int(s-win-1, x, min, max, mass, intensity, scan_index, point_count);
		
		/* compute differences */
		for (i = 0; i < N; i++)
			diff[s + i * (scan_number-1) ] = y[i] - x[i];
	}

	/* find the peak in every window */
	for (s = 1; s < scan_number - 3; s++) {
		for (i = 0; i < N; i++) {
			if ( (diff[s + i*(scan_number-1)] > 0 && diff[(s+1) + i*(scan_number-1)] < 0) ||
				 (diff[s + i*(scan_number-1)] > 0 && diff[(s+1) + i*(scan_number-1)] == 0 &&
				  diff[s+2 + i*(scan_number-1)] < 0)
				) {
				/* look for the local max within plus minus LM_WIN */
				max_int = get_max_int(s, scan_number, i+min, mass,
					intensity, scan_index, point_count, &s_max_int);

				if (max_int >= min_int)
					maxm[i * scan_number + s_max_int] = max_int;
			}
		}
	}

	UNPROTECT(6);
	return MaxIntMatrix;
}

/* Get an array of intensities in a given time point s. */

int get_int(int s, int *tmp, int min, int max, int *mass, int *intensity, 
	int * scan_index, int * point_count)
{
	int j;
	for (j = 0; j < point_count[s]; j++) {
		if(mass[scan_index[s]+j] >= min && mass[scan_index[s]+j] <= max)
			tmp[mass[scan_index[s] + j] - min] = intensity[scan_index[s] + j];
	}
	return 1;
}

/*Get the max intensity of mass "ms" and its time index in a given window (LM_WIN)*/

int get_max_int(int s, int scan_number, int ms, int *mass, int *intensity,
 int * scan_index, int *point_count, int *p)
{
	int i, j, max_int = 0, rt_max_idx = s, w = LM_WIN;
	
	for (i = -w; i <= w; i++) {
		if (s+i < 0 || s+i >= scan_number)
			continue;
			
		for (j = 0; j < point_count[s+i]; j++)
			if ( ms == mass[scan_index[s+i] + j] &&
				intensity[scan_index[s+i] + j] > max_int ) {
				
				max_int = intensity[scan_index[s+i] + j];
				rt_max_idx = s+i;
			}
	}
	*p = rt_max_idx;
	return max_int;
}
