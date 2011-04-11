/* Change Log: 
 * - 2008.06.06 Fixed: NAs were not handled correctly
*/

#include <R.h>
#include <Rdefines.h>
#include <ctype.h>
#include "getLine.h"
#include "find_peaks.h"

/* Main program */

SEXP FindPeaks(SEXP MyFile, SEXP RI_Min, SEXP Mass, SEXP RI_Max, SEXP columns) {
	FILE *fp;
	int n_scans = 0, libtotal = 0, j;
	int *mass, *int_found;
	double *ri_min, *ri_max, *ri_found, *rt_found;
	int ri_COL, sp_COL, rt_COL;
	char *myfile;
	SPECTRA *spectra;
	SEXP result, RI_Found, RT_Found, INT_Found;

	/* R objects created in the C code have to be reported using the PROTECT 
	 * macro on a pointer to the object. This tells R that the object is in 
	 * use so it is not destroyed. */
	PROTECT(MyFile  = AS_CHARACTER(MyFile));
	PROTECT(columns = AS_INTEGER(columns));
	PROTECT(RI_Min  = AS_NUMERIC(RI_Min));
	PROTECT(RI_Max  = AS_NUMERIC(RI_Max));
	PROTECT(Mass    = AS_INTEGER(Mass));

	/* Copy protected elements to pointers */
	myfile = R_alloc(strlen(CHAR(STRING_ELT(MyFile, 0))), sizeof(char));
	strcpy(myfile, CHAR(STRING_ELT(MyFile, 0)));

	ri_min = NUMERIC_POINTER(RI_Min);
	ri_max = NUMERIC_POINTER(RI_Max);
	mass   = INTEGER_POINTER(Mass);
	sp_COL = INTEGER_POINTER(columns)[0]; /* Spectra column number */
	ri_COL = INTEGER_POINTER(columns)[1]; /* R.I. column number */
	rt_COL = INTEGER_POINTER(columns)[2]; /* R.T. column number */

	libtotal = GET_LENGTH(Mass);

	/****************************************/

	fp = fopen(myfile, "r");
	is_open(fp, myfile, 5);
	spectra = read_file(fp, &n_scans, sp_COL, ri_COL, rt_COL);
	fclose(fp);

	PROTECT(RI_Found  = NEW_NUMERIC(libtotal));
	ri_found = NUMERIC_POINTER(RI_Found);
	PROTECT(RT_Found  = NEW_NUMERIC(libtotal));
	rt_found = NUMERIC_POINTER(RT_Found);
	PROTECT(INT_Found = NEW_INTEGER(libtotal));
	int_found = INTEGER_POINTER(INT_Found);

	for (j = 0; j < libtotal; j++) {
		if (ISNAN(ri_min[j]) || ISNAN(mass[j]) || ISNAN(ri_max[j])) {
			ri_found[j]  = NA_REAL;
			rt_found[j]  = NA_REAL;
			int_found[j] = NA_INTEGER;
			continue;
		}
		
		find_peak(ri_min[j], mass[j], ri_max[j], spectra, n_scans,
			ri_found+j, int_found+j, rt_found+j);
	}
	
	/* free_spectra(spectra, n_scans); */

	/* Creating a list with 2 vector elements: ri_found and int_found */
	PROTECT(result = allocVector(VECSXP, 3));
	// Attaching elements
	SET_VECTOR_ELT(result, 0, INT_Found);
	SET_VECTOR_ELT(result, 1, RI_Found);
	SET_VECTOR_ELT(result, 2, RT_Found);

	UNPROTECT(9);
	return result;
}

/* end main program */

/* Functions */

SPECTRA *read_file(FILE *fp, int *ptotal, int SPECTRUM_COL, int RI_COL, int RT_COL) {
	SPECTRA *spectra;
	int total = 0, i, j;
	int header = 1;
	char *line = NULL;
	int  len = 0;

	char *ri_str = NULL, *sp_str = NULL, *rt_str = NULL;
	int  ri_i, sp_i, rt_i, tabs, n;
	int  ri_len = 0, sp_len = 0, rt_len = 0;

	while (getLine(&line, &len, fp) != -1) {
		total++;
	}

	total--; /* header isn't counted */

	spectra = (SPECTRA *) R_alloc(total, sizeof(SPECTRA));

	fseek(fp, 0, SEEK_SET);
	j = 0;

	while (getLine(&line, &len, fp) != -1) {
		if (header) {
			header = 0;
			continue;
		}
		tabs = 0;
		i    = 0;
		ri_i = 0;
		rt_i = 0;
		sp_i = 0;
		n    = 0;

		/* allocates memory for RI string and spectra if 'line' 
		 * length is updated. string lengths will be the same as 'line' */
		str_alloc(ri_str, ri_len, len);
		str_alloc(rt_str, rt_len, len);
		str_alloc(sp_str, sp_len, len);

		while (i < strlen(line)) {
			if (line[i] == '\t' || line[i] == '\n' || line[i] == '\r')
				tabs++;
			if (tabs == RT_COL)
				if(line[i] != '\t')
					rt_str[rt_i++] = line[i];
			if (tabs == SPECTRUM_COL) {
				if(line[i] != '\t')
					sp_str[sp_i++] = line[i];
				if (line[i] == ':')
					n++;
			}
			if (tabs == RI_COL)
				if(line[i] != '\t')
					ri_str[ri_i++] = line[i];
			i++;
		}
		ri_str[ri_i] = '\0';
		rt_str[rt_i] = '\0';
		sp_str[sp_i] = '\0';
		
		if(n == 0)
			error("Error reading spectra. Invalid spectrum format\n");

		spectra[j].n  = n;
		spectra[j].ri = atof(ri_str);
		spectra[j].rt = atof(rt_str);

		spectra[j].mass = (int *) R_alloc(n , sizeof(int));
		spectra[j].in   = (int *) R_alloc(n , sizeof(int));

		if(read_spectrum(sp_str, spectra[j].mass, spectra[j].in, n) == 0)
			error("Error reading spectra. Invalid spectrum format\n");
		j++;
	}

	if(line)
		R_chk_free(line);
	if(ri_str)
		R_chk_free(ri_str);
	if(rt_str)
		R_chk_free(rt_str);
	if(sp_str)
		R_chk_free(sp_str);

	*ptotal = total;
	return spectra;
}

int read_spectrum(char *spectrum, int *mass, int *in, int n) {
	int i = 0, j = 0;
	char mass_str[STRLEN], in_str[STRLEN];
	int  mass_i = 0, in_i = 0;
	int flag = 0, sp_len;

	sp_len = strlen(spectrum);

	for (i = 0; i < sp_len; i++) {
		if (spectrum[i] == ':') {
			flag = 1;
			if(mass_i == 0) /* string is empty */
				return 0;

			mass_str[mass_i] = '\0';
			mass[j] = atoi(mass_str);
			mass_i = 0;
			continue;
		}

		if (spectrum[i] == ' ') {
			flag = 0;
			if(in_i == 0) /* string is empty */
				return 0;

			in_str[in_i] = '\0';
			in[j++] = atoi(in_str);
			in_i = 0;
			continue;
		}

		if(!isdigit(spectrum[i]))
			return 0;

		if (flag == 0)
			mass_str[mass_i++] = spectrum[i];

		if (flag == 1)
			in_str[in_i++] = spectrum[i];

		if(mass_i > STRLEN - 1 || in_i > STRLEN - 1)
			return 0;
	}
	in_str[in_i] = '\0';
	in[j] = atoi(in_str);

	return 1;
}

void find_peak
(double ri_min,int mass,double ri_max,SPECTRA *sp,int n_scans,double *ri_found, 
	int *int_found, double *rt_found) {
	int i, j, imax = -1;
	int max = -1;

/* I assume that RI is a linear function of the scan index and the RIs are
 * in ascendent order. This will speed up the search for the RI window */

	int m, b; /* slope and intercept */
	m = (int) (sp[n_scans-1].ri - sp[0].ri) / (n_scans - 1);
	b = sp[0].ri;

	/* This returns approximately the scan index of ri_min */
	i = (ri_min - b) / m;

	/* check array limits */
	i = i < 0 ? 0 : i;
	i = i > n_scans - 1 ? n_scans - 1 : i;

	while (i > 0 && ri_min < sp[i].ri)
		i--;

	for (; i < n_scans; i++) {
		if (ri_min < sp[i].ri && ri_max > sp[i].ri) {
			for (j = 0; j < sp[i].n; j++) {
				if (mass == sp[i].mass[j] && max < sp[i].in[j]) {
					max = sp[i].in[j];
					imax = i;
				}
			}
		}
		else if (ri_max < sp[i].ri)
			break;
	}
	if (imax != -1) {
		*ri_found  = sp[imax].ri;
		*rt_found  = sp[imax].rt;
		*int_found = max;
	}
	else {
		*ri_found  = NA_REAL;
		*rt_found  = NA_REAL;
		*int_found = NA_INTEGER;
	}
}

