/*
 * spectrum.c - functions to manage spectrum type
 */

#include <R.h>
#include <errno.h>
#include "spectrum.h"

#define CHK_NULL(x)  if((x) == NULL) goto error

int spectrum_add(spectrum_t *z, int mz, int in)
{
	if(z->alloc == 0) {
		z->alloc = 1024;
		z->len = 0;
		z->mz = R_Calloc(z->alloc, int);
		z->in = R_Calloc(z->alloc, int);
		if(z->mz == NULL || z->in == NULL) {
			R_Free(z->mz);
			R_Free(z->in);
			return 0;
		}
	} else if (z->len >= z->alloc) {
		int *x = R_Realloc(z->mz, 2 * z->alloc, int);
		int *y = R_Realloc(z->in, 2 * z->alloc, int);
		if(x == NULL || y == NULL) {
			R_Free(x);
			R_Free(y);
			return 0;
		}
		z->alloc *= 2;
		z->mz = x;
		z->in = y;
	}
	z->mz[z->len] = mz;
	z->in[z->len++] = in;
	return 1;
}

int spectrum_init(spectrum_t *z, int alloc)
{
	if(alloc < 0)
		return 0;
	memset(z, 0, sizeof(spectrum_t));
	if(alloc > 0) {
		CHK_NULL(z->mz = R_Calloc(alloc, int));
		CHK_NULL(z->in = R_Calloc(alloc, int));
		z->alloc = alloc;
	}
	return 1;
error:
	R_Free(z->mz);
	R_Free(z->in);
	memset(z, 0, sizeof(spectrum_t));
	return 0;
}

void _spectra_free(spectra_t * spectra)
{
	if(!spectra)
		return;
	for(int k = 0; k < spectra->n_scans; k++) {
		R_Free(spectra->sp[k].mz);
		R_Free(spectra->sp[k].in);
	}
	R_Free(spectra->sp);
	R_Free(spectra->ri);
	R_Free(spectra->rt);
	R_Free(spectra);
}

spectra_t * spectra_init(int alloc)
{
	if(alloc < 0)
		return NULL;
	spectra_t * spectra = R_Calloc(1, spectra_t);

	CHK_NULL(spectra);
	if(alloc > 0) {
		CHK_NULL(spectra->ri = R_Calloc(alloc, double));
		CHK_NULL(spectra->rt = R_Calloc(alloc, double));
		CHK_NULL(spectra->sp = R_Calloc(alloc, spectrum_t));
	}
	spectra->alloc = alloc;
	return spectra;
error:
	spectra_free(spectra);
	return NULL;
}

int spectra_realloc(spectra_t * spectra)
{
	double *ri, *rt;
	spectrum_t *sp;
	int alloc = spectra->alloc ? spectra->alloc * 2 : 1024;

	if((ri = R_Realloc(spectra->ri, alloc, double)) != NULL)
		spectra->ri = ri;
	if((rt = R_Realloc(spectra->rt, alloc, double)) != NULL)
		spectra->rt = rt;
	if((sp = R_Realloc(spectra->sp, alloc, spectrum_t)) != NULL)
		spectra->sp = sp;
	if(ri && rt && sp) {
		spectra->alloc = alloc;
		return 1;
	}
	return 0;
}

int spectra_add(spectra_t * spectra, double rt, double ri, spectrum_t *z)
{
	if(spectra->n_scans >= spectra->alloc)
		if(spectra_realloc(spectra) == 0)
			return 0;

	memcpy(spectra->sp + spectra->n_scans, z, sizeof(*z));
	spectra->ri[ spectra->n_scans ] = ri;
	spectra->rt[ spectra->n_scans ] = rt;
	spectra->n_scans++;
	return 1;
}

#define VALINT(x) ( (x) >= 0 && (x) <= 0x7FFFFFFF )
#define NULLOREMPTY(x) ( (x == NULL) || ( *x == '\0' ) )

#include "spectrum_re.h" /* generated code from spectrum.re */

/**
 * Parses a spectrum string and returns the number of m/z - int pairs
 *
 * The spectrum string is a list of m/z-intensity values separated by spaces.
 * The m/z and intensity values must be integers and be separated by a single
 * colon (:). Separation by multiple spaces, as well as leading and trailing
 * spaces are ignored.
 *
 * The function expects a pointer `spectrum_t` struct which is used to store
 * the values.
 *
 * @param str  The NIL-terminated string to check validity
 * @param z    A pointer to the `spectrum_t` object
 * @return On success, the number of m/z - intensity pairs reads (a positive
 *     value). On failure, it returns a negative value:
 *
 *     -1: An invalid character was detected.
 *     -2: Unable to convert a string to integer.
 *     -3: A error was raised by sscanf(3) or and integer overflow was detected.
 *     -4: Unable to allocate memory in the spectrum object
 *
 *     The special value of zero (0) means that the string was empty or contained
 *     only whitespaces.
 *
 * @note It is possible that some m/z-integer pairs are assigned despite a
 *      failure. This is because the values are saved as soon as their are read.
 *
 * generate C code:
 * re2c <input> -o <output> -i --tags
 *
 * reference: https://re2c.org/manual/manual_c.html#sentinel
 */

int spectrum_parse(const char *str, spectrum_t * z)
{
	if(NULLOREMPTY(str))
		return 0;

	const char *start = NULL, *end = NULL;
	int len = z->len, count = 0, ret = 0;

	while((ret = _spectrum_re(str, &start, &end)) != 1) {
		long m, i;
		if(ret == SP_ERR_INVALID)
			goto error;
		errno = 0;
		if(sscanf(start, "%ld:%ld", &m, &i) != 2) {
			ret = SP_ERR_STRTOL;
			goto error;
		}

		if(errno != 0 || !VALINT(m) || !VALINT(i)) {
			ret = SP_ERR_OVERFLOW;
			goto error;
		}
		if(spectrum_add(z, (int) m, (int) i) == 0) {
			ret = SP_ERR_ALLOC;
			goto error;
		}
		str = end;
		count++;
	}
	return count;
error:
	z->len = len;
	return ret;
}

/* scan string and return the number of ':' of potential separators */
int spectrum_scan(const char *str)
{
	int count = 0;
	while(*str) {
		if(*str == ':')
			count++;
		str++;
	}
	return count;
}
