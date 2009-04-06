/* Change Log:
	- 2008.06.09: Fixed. MassLimits were not properly checked.
*/

#include <R.h>
#include <Rdefines.h>

/* Extracts all CDF Intensities into a Matrix */

/*  Main Program */

SEXP peakExtraction(SEXP MassValues, SEXP IntensityValues, SEXP PointCount,
	SEXP ScanIndex, SEXP MassLimits) {
	
	int min, max, N;
	int scan_number, *scan_index, *point_count, *intensity, *mass;
	
	int s, i, *maxm;
	
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
	
	N = max - min + 1;
	
	scan_number = GET_LENGTH(ScanIndex);
	
	scan_index  = INTEGER_POINTER(ScanIndex);
	point_count = INTEGER_POINTER(PointCount);
	intensity   = INTEGER_POINTER(IntensityValues);
	mass        = INTEGER_POINTER(MassValues);

	PROTECT(MaxIntMatrix = allocMatrix(INTSXP, scan_number, N));
	maxm = INTEGER_POINTER(MaxIntMatrix);
	
	for(s = 0; s < scan_number; s++) {
		 /* Make sure that the matrix is filled with zeros */
                for (i = 0; i < N; i++)
                        maxm[i*scan_number+s] = 0;

		for (i = 0; i < point_count[s]; i++)
			if(mass[scan_index[s]+i] >= min && mass[scan_index[s]+i] <= max)
				maxm[(mass[scan_index[s] + i]-min) * scan_number + s] = intensity[scan_index[s] + i];

	}
	
	UNPROTECT(6);
	return MaxIntMatrix;
}

