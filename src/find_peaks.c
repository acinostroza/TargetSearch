/* find peaks in dat files */
#include <R.h>
#include <Rdefines.h>
#include <ctype.h>

#include "find_peaks.h"
#include "getLine.h"
#include "utils.h"

#define BUFFER 4096

/* search type: */
typedef enum {
	ALL=1,   /* return all peaks */
	MINRI,   /* return peaks with min RI error */
	MAXINT   /* return most abundant peaks */
} SearchType;

struct point_list_s *
init_point_list(int size)
{
	struct point_list_s *x = Calloc(1, struct point_list_s);
	x->p = Calloc(size, struct point_s);
	x->alloc = size;
	x->length = 0;
	return x;
}

void
add_point(struct point_list_s *x, struct point_s *p)
{
	int k = x->length;
	if(k == x->alloc) {
		x->alloc *= 2;
		x->p = Realloc(x->p, x->alloc, struct point_s);
	}
	x->p[k] = *p;
	x->length++;
}

void
free_point_list(struct point_list_s *x)
{
	Free(x->p);
	Free(x);
}

void
find_all_peaks(double mass, double ri_exp, double ri_min,  double ri_max,
		SPECTRA *sp, struct point_list_s *plist, int use_rt, int idx)
{
	int i, j;
	int n_scans = sp->n_scans;
	double *ri;
	struct point_s p;

	/* uses RT or RI to perform the search */
	ri = (use_rt == 0) ? sp->ri : sp->rt;

	/* This returns approximately the scan index of ri_min */
	i = binsearch(ri, ri_min, n_scans);

	for (; i < n_scans; i++) {
		if (ri_max < ri[i])
			break;
		else if (ri_min < ri[i] && ri_max > ri[i]) {
			for (j = 0; j < sp->n[i]; j++) {
				if (mass == sp->pk[i].mass[j]) {
					p.rt = sp->rt[i];
					p.ri = sp->ri[i];
					p.in = sp->pk[i].in[j];
					p.mz = sp->pk[i].mass[j];
					p.idx = idx;
					p.err = fabs(ri_exp - ri[i]);
					add_point(plist, &p);
				}
			}
		}
	}
}

struct point_list_s *
filter_results(struct point_list_s *plist, SearchType st)
{
	struct point_list_s *result;
	int prev = -1;
	struct point_s *best = NULL;

	if(st == ALL || plist->length <= 1)
		return plist;

	result = init_point_list(plist->length);

	for(int i = 0; i < plist->length; i++) {
		struct point_s *p = plist->p + i;
		if(p->idx != prev) {
			if(best != NULL) {
				add_point(result, best);
			}
			prev = p->idx;
			best = p;
			continue;
		}
		if(st == MAXINT) {
			if (best->in < p->in)
				best = p;
			continue;
		}
		if(st == MINRI) {
			if (best->err > p->err)
				best = p;
			continue;
		}
	}
	if(best != NULL) {
		add_point(result, best);
	}
	return result;
}

SPECTRA *
read_file(const char *file, int ftype, int swap, int sp_COL, int ri_COL, int rt_COL)
{
	FILE *fp;
	SPECTRA *spectra;
	if(ftype == 0) {
		fp = fopen(file, "r");
		if(fp == NULL)
			error("Error opening file %s\n", file);
		spectra = read_txt(fp, sp_COL, ri_COL, rt_COL);
		if(!spectra)
			error("Error reading file %s\n", file);
	} else {
		fp = fopen(file, "rb");
		if(fp == NULL)
			error("Error opening file %s\n", file);
		spectra = read_dat(fp, swap);
		if(!spectra)
			error("Error reading file %s\n", file);
	}
	fclose(fp);
	return spectra;
}

struct point_list_s *
do_search(SPECTRA *spectra, int *mass, double *ri_exp, double *ri_min, double *ri_max,
		int use_rt, SearchType st, int libtotal)
{
	struct point_list_s *plist = init_point_list(2 * libtotal);
	struct point_list_s *res;

	for (int j = 0; j < libtotal; j++) {
		double ri = (ri_exp == NULL) ? 0.0 : ri_exp[j];
		if (ISNAN(ri_min[j]) || (mass[j] == NA_INTEGER) || ISNAN(ri_max[j]))
			continue;
		find_all_peaks(mass[j], ri, ri_min[j], ri_max[j], spectra, plist, use_rt, j);
	}
	res = filter_results(plist, st);
	if(plist != res)
		free_point_list(plist);
	return res;
}

/*
  Find Peaks version 2.
  Description:
     Looks for peaks in a RI file. It returns either (1) all peaks found within
     the time range, (2) the closest to the expected RI, or (3) the most intense.
  Output:
     1. vector of intensities.
     2. vector of retention indexes.
     3. vector of retention times
     4. vector of indexes relative to the input.
  Args:
     RI_file. RI file (dat or txt)
     Mass. vector of m/z values to search.
     RIexp. Expected Retention Index (RI)
     RI_Min. RI minimum value.
     RI_Max. RI maximum value.
     Options. Vector of file format options (see comment below).
     useRT. Should use RT or RI. Obviouly RIexp, RI_Min, RI_max units should be
            Time in seconds or Index units.
     Search: 1. all peaks, 2. closest to expected RI, 3. most intense.
  Notes:
     RI_exp can be NULL, but filterResult cannot be 1. In this case RI_exp
     is ignored.
 */

SEXP find_peaks(SEXP RI_file, SEXP Mass, SEXP RI_exp, SEXP RI_Min, SEXP RI_Max, SEXP Options,
	SEXP useRT, SEXP Search)
{
	/* internal variables */
	char *file;
	int  libtotal = GET_LENGTH(Mass);
	SPECTRA *spectra;
	struct point_list_s *res;

	/* R variables */
	int *mass  = INTEGER(Mass);
	int use_rt = LOGICAL(useRT)[0];
	SearchType st  = (SearchType) INTEGER(Search)[0];
	double *ri_min = REAL(RI_Min), *ri_max = REAL(RI_Max);
	double *ri_exp = isNull(RI_exp) ? NULL : REAL(RI_exp);

	/* output variables */
	SEXP RI_Found, RT_Found, INT_Found, I_Found, result;

	/* parse options */
	int ftype  = INTEGER(Options)[0]; /* file type: 0 = TXT; 1 = DAT */
	int swap   = INTEGER(Options)[1]; /* swap = 1 in big endian platforms */
	int sp_COL = INTEGER(Options)[2]; /* Spectra column number */
	int ri_COL = INTEGER(Options)[3]; /* R.I. column number */
	int rt_COL = INTEGER(Options)[4]; /* R.T. column number */

	/* Copy file name to a string */
	file = R_alloc(strlen(CHAR(STRING_ELT(RI_file, 0))), sizeof(char));
	strcpy(file, CHAR(STRING_ELT(RI_file, 0)));

	/* parse file */
	spectra = read_file(file, ftype, swap, sp_COL, ri_COL, rt_COL);

	/* search peaks */
	res = do_search(spectra, mass, ri_exp, ri_min, ri_max, use_rt, st, libtotal);

	RI_Found  = PROTECT(allocVector(REALSXP, res->length));
	RT_Found  = PROTECT(allocVector(REALSXP, res->length));
	INT_Found = PROTECT(allocVector(INTSXP, res->length));
	I_Found   = PROTECT(allocVector(INTSXP, res->length));

	for(int j = 0; j < res->length; j++) {
		struct point_s *p = res->p + j;
		REAL(RI_Found)[j] = p->ri;
		REAL(RT_Found)[j] = p->rt;
		INTEGER(INT_Found)[j] = p->in;
		INTEGER(I_Found)[j] = p->idx;
	}
	/* Creating a list with 2 vector elements: ri_found and int_found */
	PROTECT(result = allocVector(VECSXP, 4));
	// Attaching elements
	SET_VECTOR_ELT(result, 0, INT_Found);
	SET_VECTOR_ELT(result, 1, RI_Found);
	SET_VECTOR_ELT(result, 2, RT_Found);
	SET_VECTOR_ELT(result, 3, I_Found);

	free_point_list(res);
	UNPROTECT(5);
	return result;
}

