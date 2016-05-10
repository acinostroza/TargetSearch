/* find peaks in dat files */
#include <R.h>
#include <Rdefines.h>
#include <ctype.h>

#include "file.h"
#include "find_peaks.h"
#include "getLine.h"

/* Main program */

SEXP FindPeaks(SEXP MyFile, SEXP RI_Min, SEXP Mass, SEXP RI_Max, SEXP Options,
	SEXP useRT)
{
        FILE *fp;
        int  libtotal = 0, j;
        int *mass, *int_found;
        double *ri_min, *ri_max, *ri_found, *rt_found;
        char *myfile;
        SPECTRA *spectra;
        SEXP result, RI_Found, RT_Found, INT_Found;

	int ftype; /* file type: 0 = TXT; 1 = DAT */
	int ri_COL, sp_COL, rt_COL;
	int swap;
	int use_rt; /* use RT instead of RI */

        /* R objects created in the C code have to be reported using the PROTECT
         * macro on a pointer to the object. This tells R that the object is in
         * use so it is not destroyed. */
        PROTECT(MyFile  = AS_CHARACTER(MyFile));
        PROTECT(RI_Min  = AS_NUMERIC(RI_Min));
        PROTECT(RI_Max  = AS_NUMERIC(RI_Max));
        PROTECT(Mass    = AS_INTEGER(Mass));
        PROTECT(Options = AS_INTEGER(Options));
        PROTECT(useRT   = AS_INTEGER(useRT));

        /* Copy protected elements to pointers */
        myfile = R_alloc(strlen(CHAR(STRING_ELT(MyFile, 0))), sizeof(char));
        strcpy(myfile, CHAR(STRING_ELT(MyFile, 0)));

        ri_min = NUMERIC_POINTER(RI_Min);
        ri_max = NUMERIC_POINTER(RI_Max);
        mass   = INTEGER_POINTER(Mass);

	ftype  = INTEGER_POINTER(Options)[0]; /* file type: 0 = TXT; 1 = DAT */
	swap   = INTEGER_POINTER(Options)[1]; /* swap = 1 in big endian platforms */
        sp_COL = INTEGER_POINTER(Options)[2]; /* Spectra column number */
        ri_COL = INTEGER_POINTER(Options)[3]; /* R.I. column number */
        rt_COL = INTEGER_POINTER(Options)[4]; /* R.T. column number */
	use_rt = INTEGER_POINTER(useRT)[0];   /* use RT instead of RI */

        libtotal = GET_LENGTH(Mass);

        /****************************************/

	if(ftype == 0) {
		fp = fopen(myfile, "r");
		is_open(fp, myfile, 5);
		spectra = read_txt(fp, sp_COL, ri_COL, rt_COL);
		if(!spectra)
			error("Error reading file %s\n", myfile);
	} else {
		fp = fopen(myfile, "rb");
		is_open(fp, myfile, 5);
		spectra = read_dat(fp, swap);
		if(!spectra)
			error("Error reading file %s\n", myfile);
	}
        fclose(fp);

        PROTECT(RI_Found  = NEW_NUMERIC(libtotal));
        ri_found = NUMERIC_POINTER(RI_Found);
        PROTECT(RT_Found  = NEW_NUMERIC(libtotal));
        rt_found = NUMERIC_POINTER(RT_Found);
        PROTECT(INT_Found = NEW_INTEGER(libtotal));
        int_found = INTEGER_POINTER(INT_Found);

        for (j = 0; j < libtotal; j++) {
                if (ISNAN(ri_min[j]) || (mass[j] == NA_INTEGER) || ISNAN(ri_max[j])) {
                        ri_found[j]  = NA_REAL;
                        rt_found[j]  = NA_REAL;
                        int_found[j] = NA_INTEGER;
                        continue;
                }

                find_peak(ri_min[j], mass[j], ri_max[j], spectra,
                        ri_found+j, int_found+j, rt_found+j, use_rt);
        }

        /* free_spectra(spectra, n_scans); */

        /* Creating a list with 2 vector elements: ri_found and int_found */
        PROTECT(result = allocVector(VECSXP, 3));
        // Attaching elements
        SET_VECTOR_ELT(result, 0, INT_Found);
        SET_VECTOR_ELT(result, 1, RI_Found);
        SET_VECTOR_ELT(result, 2, RT_Found);

        UNPROTECT(10);
        return result;
}

void find_peak
(double ri_min,int mass,double ri_max,SPECTRA *sp, double *ri_found,
        int *int_found, double *rt_found, int use_rt)
{
        int i, j, imax = -1;
        int max = -1;
	int n_scans = sp->n_scans;
	double *ri;

/* I assume that RI is a linear function of the scan index and the RIs are
 * in ascendent order. This will speed up the search for the RI window */

	/* uses RT or RI to perform the search */
	ri = (use_rt == 0) ? sp->ri : sp->rt;

        int m, b; /* slope and intercept */
        m = (int) (ri[n_scans-1] - ri[0]) / (n_scans - 1);
        b = ri[0];

        /* This returns approximately the scan index of ri_min */
        i = (ri_min - b) / m;

        /* check array limits */
        i = i < 0 ? 0 : i;
        i = i > n_scans - 1 ? n_scans - 1 : i;

        while (i > 0 && ri_min < ri[i])
                i--;

        for (; i < n_scans; i++) {
                if (ri_min < ri[i] && ri_max > ri[i]) {
                        for (j = 0; j < sp->n[i]; j++) {
                                if (mass == sp->pk[i].mass[j] && max < sp->pk[i].in[j]) {
                                        max = sp->pk[i].in[j];
                                        imax = i;
                                }
                        }
                }
                else if (ri_max < ri[i])
                        break;
        }
        if (imax != -1) {
                *ri_found  = sp->ri[imax];
                *rt_found  = sp->rt[imax];
                *int_found = max;
        }
        else {
                *ri_found  = NA_REAL;
                *rt_found  = NA_REAL;
                *int_found = NA_INTEGER;
        }
}

