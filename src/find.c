/* find peaks in dat files */
#include <R.h>
#include <Rdefines.h>
#include <ctype.h>

#include "find.h"
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
	struct point_list_s *x = R_Calloc(1, struct point_list_s);
	if(x == NULL)
		return NULL;
	if((x->p = R_Calloc(size, struct point_s)) == NULL) {
		R_Free(x);
		return NULL;
	}
	x->alloc = size;
	x->length = 0;
	return x;
}

int
add_point(struct point_list_s *x, struct point_s *p)
{
	int k = x->length;
	if(k == x->alloc) {
		int alloc = x->alloc ? 2 * x->alloc : BUFFER;
		struct point_s * q = R_Realloc(x->p, alloc, struct point_s);
		if(q == NULL)
			return 0;
		x->p = q;
		x->alloc = alloc;
	}
	x->p[k] = *p;
	x->length++;
	return 1;
}

void
free_point_list(struct point_list_s *x)
{
	if(x) {
		R_Free(x->p);
		R_Free(x);
	}
}

int
find_all_peaks(double mass, double ri_exp, double ri_min,  double ri_max,
		spectra_t *sp, struct point_list_s *plist, int use_rt, int idx)
{
	int i, j;
	struct point_s p;

	/* uses RT or RI to perform the search */
	double *ri = (use_rt == 0) ? sp->ri : sp->rt;

	/* This returns approximately the scan index of ri_min */
	i = binsearch(ri, ri_min, sp->n_scans);

	for (; i < sp->n_scans; i++) {
		if (ri_max < ri[i])
			break;
		else if (ri_min < ri[i] && ri_max > ri[i]) {
			for (j = 0; j < sp->sp[i].len; j++) {
				if (mass == sp->sp[i].mz[j]) {
					p.rt = sp->rt[i];
					p.ri = sp->ri[i];
					p.in = sp->sp[i].in[j];
					p.mz = sp->sp[i].mz[j];
					p.idx = idx;
					p.err = fabs(ri_exp - ri[i]);
					if(!add_point(plist, &p))
						return 0;
				}
			}
		}
	}
	return 1;
}

struct point_list_s *
filter_results(struct point_list_s *plist, SearchType st)
{
	struct point_list_s *result;
	int prev = -1;
	struct point_s *best = NULL;

	if(st == ALL || plist->length <= 1)
		return plist;

	if((result = init_point_list(plist->length)) == NULL)
		return NULL;

	for(int i = 0; i < plist->length; i++) {
		struct point_s *p = plist->p + i;
		if(p->idx != prev) {
			if(best != NULL) {
				if(!add_point(result, best))
					goto error;
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
		if(!add_point(result, best))
			goto error;
	}
	return result;
error:
	free_point_list(result);
	return NULL;
}

spectra_t *
read_file(const char *file, int ftype, int swap,
	  const char *sp_COL, const char *ri_COL, const char *rt_COL, const int * cols)
{
	FILE *fp = NULL;
	char *mode = ftype == 0 ? "rt" : "rb";

	if((fp = fopen(file, mode)) == NULL)
		error("Error opening file %s\n", file);

	spectra_t * spectra = ftype == 0 ? read_txt(fp, sp_COL, ri_COL, rt_COL, cols) :
				read_dat(fp, swap);
	fclose(fp);
	if(!spectra)
		error("Error reading file %s\n", file);
	return spectra;
}

struct point_list_s *
do_search(spectra_t *spectra, int *mass, double *ri_exp, double *ri_min, double *ri_max,
		int use_rt, SearchType st, int libtotal)
{
	struct point_list_s *plist = init_point_list(2 * libtotal);
	struct point_list_s *res = NULL;
	if(plist == NULL)
		return NULL;

	for (int j = 0; j < libtotal; j++) {
		double ri = (ri_exp == NULL) ? 0.0 : ri_exp[j];
		if (ISNAN(ri_min[j]) || (mass[j] == NA_INTEGER) || ISNAN(ri_max[j]))
			continue;
		if(!find_all_peaks(mass[j], ri, ri_min[j], ri_max[j], spectra, plist, use_rt, j))
			goto error;
	}
	if((res = filter_results(plist, st)) == NULL)
		goto error;
	if(plist != res)
		free_point_list(plist);
	return res;
error:
	free_point_list(plist);
	free_point_list(res);
	return NULL;
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
	const char *file = CHARACTER_VALUE(RI_file);
	int  libtotal = GET_LENGTH(Mass);
	spectra_t *spectra;
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
	int cols[] = {sp_COL, ri_COL, rt_COL};

	/* parse file */
	spectra = read_file(file, ftype, swap, NULL, NULL, NULL, cols);

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
