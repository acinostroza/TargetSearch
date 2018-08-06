#ifndef _FIND_PEAKS_H
#define _FIND_PEAKS_H

/* require SPECTRA data type */
#include "file.h"

/* Function prototypes */

struct point_s {
	double rt;   /* Ret. Time  */
	double ri;   /* Ret. Index */
	int    in;   /* Intensity  */
	int    mz;   /* mass to charge */
	double err;  /* abs time error */
	int    idx;  /* library index */
};

struct point_list_s {
	struct point_s *p;
	int length;
	int alloc;
};

void
find_all_peaks(double, double, double, double, SPECTRA *,
		struct point_list_s *, int, int);

SEXP find_peaks(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#endif
