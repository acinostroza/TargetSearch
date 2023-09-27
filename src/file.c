/* Functions to perform file manipulation */
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <ctype.h>
#include "file.h"
#include "strutils.h"
#include "get_line.h"

/* a random signature for the dat files */
#define SIGLEN 8  /* signature length */

const unsigned char signature[]={42,243,27,10,183,75,1,0};

/* from R source code: swap bytes for big endian machines */
static void swapb(void *ptr, int size)
{
	int i;
	char *p=ptr, tmp;
	if(size==1)
		return;
	for(i=0; i < size/2; i++) {
		tmp  = p[i];
		p[i] = p[size - i - 1];
		p[size - i - 1] = tmp;
	}
}

/* apply byte swapping to an array */
static void swapp(void *ptr, int size, int len) {
	int i;
	char *p=ptr;
	for(i = 0; i < len; i++)
		swapb(p + i*size, size);
}

/* function to read a DAT file:
 * args:
 *   - fp: file pointer
 *   - p_splen: spectra legth
 *   - swap: 1 (big endian; swap) 0 (little endian; no swap)
 * returns: SPECTRA struct */
spectra_t * read_dat(FILE *fp, int swap)
{
	int splen=0, sptotal=0, *pcount = NULL, count = 0;
	unsigned char s[SIGLEN];
	spectra_t * spectra = NULL;

	if(fread(s, 1, SIGLEN, fp) != SIGLEN) {
		REprintf("Incorrect file signature\n");
		return NULL;
	}
	if(memcmp(s, signature, SIGLEN)) {
		REprintf("Incorrect file signature\n");
		return NULL;
	}

	if(!fread(&splen, sizeof(splen), 1, fp))
		return NULL;
	if(!fread(&sptotal, sizeof(sptotal), 1, fp))
		return NULL;

	if(swap == 1) {
		swapb(&splen, sizeof(splen));
		swapb(&sptotal, sizeof(sptotal));
	}
	if((spectra = spectra_init(splen)) == NULL)
		goto error;
	if(splen <= 0 || (pcount = R_Calloc(splen, int)) == NULL)
		goto error;

	if(!fread(spectra->ri, sizeof(double), splen, fp))
		goto error;
	if(!fread(spectra->rt, sizeof(double), splen, fp))
		goto error;
	if(!fread(pcount, sizeof(int), splen, fp))
		goto error;

	if(swap == 1) {
		swapp(pcount, sizeof(int), splen);
		swapp(spectra->ri, sizeof(double), splen);
		swapp(spectra->rt, sizeof(double), splen);
	}
	spectra->n_scans = splen;
	for(int i = 0; i < splen; i++) {
		spectrum_t * z = spectra->sp + i;
		if(!spectrum_init(z, pcount[i]))
			goto error;
		if((z->len = pcount[i]) == 0)
			continue;
		if(!fread(z->mz, sizeof(int), z->len, fp))
			goto error;
		if(!fread(z->in, sizeof(int), z->len, fp))
			goto error;
		if(swap == 1) {
			swapp(z->mz, sizeof(int), z->len);
			swapp(z->in, sizeof(int), z->len);
		}
		count += pcount[i];
	}
	if(count != sptotal)
		goto error;

	R_Free(pcount);
	return spectra;
error:
	R_Free(pcount);
	spectra_free(spectra);
	return NULL;
}

/* function to write a dat file */
int write_dat(FILE *fp, spectra_t *sp, int swap)
{
#define f_write(x) do { \
	if(fwrite(&(x), sizeof(x), 1, fp) != 1) return 0; } while(0)

        int i, j, tmp, n, splen=sp->n_scans, pcount = 0;
	double ri, rt;

	/* write files signatures */
	f_write(signature);

        /* write len*/
	n = sp->n_scans;
	if(swap == 1)
		swapb(&n, sizeof(n));
	f_write(n);

	/* write point count */
	for(i = 0; i < splen; i++)
		pcount += sp->sp[i].len;
	if(swap == 1)
		swapb(&pcount, sizeof(n));

	f_write(pcount);

        /* write RI */
        for(i = 0; i < splen; i++) {
		ri = sp->ri[i];
		if(swap == 1)
			swapb(&ri, sizeof(ri));
                f_write(ri);
	}

        /* write RT */
        for(i = 0; i < splen; i++) {
		rt = sp->rt[i];
		if(swap == 1)
			swapb(&rt, sizeof(rt));
		f_write(rt);
	}

        /* write N */
        for(i = 0; i < splen; i++) {
		n = sp->sp[i].len;
		if(swap == 1)
			swapb(&n, sizeof(n));
		f_write(n);
	}

        for(i = 0; i < splen; i++) {
		for(j = 0; j < sp->sp[i].len; j++) {
			tmp = sp->sp[i].mz[j];
			if(swap == 1)
				swapb(&tmp, sizeof(tmp));
			f_write(tmp);
		}
		for(j = 0; j < sp->sp[i].len; j++) {
			tmp = sp->sp[i].in[j];
			if(swap == 1)
				swapb(&tmp, sizeof(tmp));
			f_write(tmp);
		}
        }
	return 1;
#undef f_write
}

/* macros needed to read a TXT file with peak data */
#define MAX(a, b)  ((a) > (b) ? (a) : (b))
#define CHKCOL(x, col) do \
	if((x) < 0) { \
		REprintf("Unable to find colum `%s'\n", col ? col : "NULL"); \
		goto clean; \
	} while(0)

/**
 * Read peak data from a TXT file
 *
 * The functions reads a tab delimited file that contains peak data. It is not
 * possible to change the delimiter character as it is hard-coded.
 *
 * @param fp pointer to stream to be parsed
 * @param SPECTRUM the column name that contains the SPECTRUM data.
 * @param RI the column name that contains the retention index values.
 * @param RT the column name that contains the retention time values
 * @param cols a pointer to the numeric positions of the spectrum, RI, and RT
 *        or NULL. If not NULL, then it's size must be equal to three (this is not
 *        checked). Also in this case, the values of SPECTRUM, RI and RT are ignored.
 * @returns
 *        A spectra_t object with peak data or NULL in case of error.
 */

spectra_t * read_txt(FILE *fp, const char * SPECTRUM, const char * RI, const char * RT, const int *cols)
{
	int lncount = 0, len = 0, nread = 0, nextchr = EMPTY_CHAR, ret = 0;
	int SPECTRUM_COL = 0, RI_COL = 0, RT_COL = 0, mintab = 0;
	char * line = NULL, * temp = NULL;
	spectra_t * spectra = spectra_init(BUFSIZ);

	if(spectra == NULL) {
		REprintf("Unable to (re)allocate memory\n");
		goto clean;
	}

	/* column positions are passed explicitely */
	if(cols != NULL) {
		SPECTRUM_COL = cols[0]; RI_COL = cols[1]; RT_COL = cols[2];
		mintab = MAX(SPECTRUM_COL, MAX(RI_COL, RT_COL)) + 1;
	}

	while((nread = get_line(&line, &len, &nextchr, fp)) > 0) {
		int tabs = 0, pc = 0;
		double rt = 0, ri = 0;
		nread = rstrip(line);
		spectrum_t p = { NULL, NULL, 0, 0};

		if(++lncount == 1) {
			if(!ascii(line, (size_t) nread)) {
				REprintf("Non-ascii characters detected in header.\n");
				goto clean;
			}
			if(cols != NULL) /* do not check columns if already given */
				continue;
			CHKCOL(SPECTRUM_COL = get_col_index(line, SPECTRUM, '\t'), SPECTRUM);
			CHKCOL(RI_COL = get_col_index(line, RI, '\t'), RI);
			CHKCOL(RT_COL = get_col_index(line, RT, '\t'), RT);
			mintab = MAX(SPECTRUM_COL, MAX(RI_COL, RT_COL)) + 1;
			continue;
		}

		temp = line;
		while(temp) {
			char * next = tokenize(temp, '\t');
			if(tabs == RT_COL)
				if(!stod(temp, &rt))
					goto error;
			if(tabs == RI_COL)
				if(!stod(temp, &ri))
					goto error;
			if(tabs == SPECTRUM_COL) {
				pc = spectrum_scan(temp);
				if(!spectrum_init(&p, pc))
					goto error;
				if((ret = spectrum_parse(temp, &p)) < 0) {
					R_Free(p.in);
					R_Free(p.mz);
					goto error;
				}
			}
			temp = next;
			tabs++;
		}

		if(tabs < mintab) {
			untokenize(line, nread, '\t');
			REprintf("Not enough columns at line %d: found = %d, expected = %d\n", lncount, tabs, mintab);
			REprintf("line: %s\n", line);
			goto clean;
		}
		spectra_add(spectra, rt, ri, &p);
	}
	if(nread == ALLOC_ERROR) {
		REprintf("An error ocurred allocating memory\n");
		goto clean;
	}
	R_Free(line);
	return spectra;

error:
	REprintf("Unable to parse field at line %d: `%s`\n", lncount, temp);

clean:
	R_Free(line);
	spectra_free(spectra);
	return NULL;
}

/* transform a peak matrix to SPECTRA struct */
spectra_t * pktosp(double *rt, double *ri, int *in, int *mass, int nscans)
{
	spectra_t *sp = NULL;
	int nmz = mass[1] - mass[0] + 1;

	if(nmz < 0)
		return NULL;

	int alloc = (nmz / 5) < 32 ? 32 : (nmz / 5);

	if((sp = spectra_init(nscans)) == NULL)
		return NULL;

	sp->n_scans = nscans;

	/* construct the structure */
	for(int i = 0; i < nscans; i++) {
		spectrum_t * z = sp->sp + i;
		sp->rt[i] = rt[i];
		sp->ri[i] = ri[i];

		if(spectrum_init(z, alloc) == 0)
			goto error;

		for(int j = 0; j < nmz; j++) {
			int m = mass[0] + j, x = in[i + j * nscans];
			if(x <= 0)
				continue;
			if(spectrum_add(z, m, x) == 0)
				goto error;
		}
	}
	return sp;
error:
	spectra_free(sp);
	return NULL;
}

#define ERROR(...) do { REprintf(__VA_ARGS__); goto error; } while(0)

/* C wrapper to write_dat / write_txt functions */
static int _write_peaks(const char *fout, double *rt, double *ri, int *in, int *mass,
	int nscans, int swap, const char *header, char ftype)
{
	int ret = 0;
	spectra_t *sp = NULL;
	FILE *fp = NULL;

	if((sp = pktosp(rt, ri, in, mass, nscans)) == NULL)
		ERROR("Error creating spectra struct.\n");

	if((fp = fopen(fout, "wb")) == NULL)
		ERROR("Error writing file %s\n", fout);

	if(ftype == 'b')
		write_dat(fp, sp, swap);
	else if(ftype == 't')
		write_txt(fp, sp, header);
	ret = 1;
error:
	spectra_free(sp);
	if(fp) fclose(fp);
	return ret;
}
#undef ERROR

/* R interface to C wrapper */
SEXP write_peaks(SEXP Output, SEXP RT, SEXP RI, SEXP IN, SEXP Mass, SEXP Swap, SEXP Header)
{
	const char *file = CHARACTER_VALUE(Output), *header = NULL;
	int swap = INTEGER_VALUE(Swap), *in = INTEGER(IN), *mass = INTEGER(Mass);
	double *rt = REAL(RT), *ri = REAL(RI);
	int nscans = GET_LENGTH(RI);
	char ftype = 'b';
	SEXP ans = PROTECT(NEW_LOGICAL(1));

	if(swap == -1) {
		header = CHARACTER_VALUE(Header);
		ftype = 't';
	}
	int ret = _write_peaks(file, rt, ri, in, mass, nscans, swap, header, ftype);
	SET_LOGICAL_ELT(ans, 0, ret);
	UNPROTECT(1);
	return ans;
}

/**
 * Writes peak data to a tab-delimited file
 *
 * @param fp stream pointer to the destination
 * @param sp pointer to the peak data structure
 * @param header string containing a tab-delimited header (the convention is to
 *        write the RT column follow by the SPECTRUM and the RI columns).
 *        Example string "RT\tSPECTRUM\tRI" (note that the columns are joined
 *        and separated by a tab char).
 * @return 1 on success, 0 otherwise
 */
int write_txt(FILE *fp, spectra_t *sp, const char *header)
{
#define f_printf(...) if(fprintf(__VA_ARGS__) < 0) return 0
	f_printf(fp, "%s\n", header);
	for(int i = 0; i < sp->n_scans; i++) {
		spectrum_t *p = sp->sp + i;
		if(p->len == 0)
			continue;
		f_printf(fp, "%.15g\t", sp->rt[i]);
		for(int j = 0; j < p->len; j++)
			f_printf(fp, "%s%d:%d", j > 0 ? " " : "", p->mz[j], p->in[j]);
		f_printf(fp, "\t%.15g\n", sp->ri[i]);
	}
	return 1;
#undef f_printf
}

/* get column names or positions depending on the type of columns */
int get_columns(SEXP columns, struct column_s * col)
{
	memset(col, 0, sizeof(*col));
	if(isNull(columns) || GET_LENGTH(columns) != 3)
		return 0;
	if(isString(columns)) {
		col->sp_col = CHAR(STRING_ELT(columns, 0));
		col->ri_col = CHAR(STRING_ELT(columns, 1));
		col->rt_col = CHAR(STRING_ELT(columns, 2));
		return 1;
	}
	if(isInteger(columns)) {
		col->icols = INTEGER(columns);
		return 1;
	}
	return 0;
}

static inline SEXP mkInt(int x)
{
	SEXP ans = PROTECT(NEW_INTEGER(1));
	SET_LOGICAL_ELT(ans, 0, x);
	UNPROTECT(1);
	return ans;
}

static inline spectra_t *
rdfile(FILE *fp, int type, int swap, const struct column_s * c)
{
	return type == 0 ? read_dat(fp, swap) :
		read_txt(fp, c->sp_col, c->ri_col, c->rt_col, c->icols);
}

static inline int
wrfile(FILE *fp, spectra_t * s, int type, int swap, const char * header)
{
	return type == 0 ? write_dat(fp, s, swap) : write_txt(fp, s, header);
}

#define err(x, ...) do { ret = x; REprintf(__VA_ARGS__); goto error; } while(0)
#define close(x) if(x) fclose(x);

SEXP convert_ri_file(SEXP IN, SEXP OUT, SEXP Type, SEXP Columns, SEXP Header)
{
	FILE *fin = NULL, *fout = NULL;
	int type = INTEGER_VALUE(Type);
	const char *header = CHARACTER_VALUE(Header);
	struct column_s c;
	int swap = endianness();
	int ret = 0;
	spectra_t * sp = NULL;

	if(type == 1 && get_columns(Columns, &c) == 0)
		err(1, "Unable to parse file columns\n");
	if(!(fin = fopen(CHARACTER_VALUE(IN), type == 1 ? "rt" : "rb")))
		err(2, "Unable to open file `%s`\n", CHARACTER_VALUE(IN));
	if(!(fout = fopen(CHARACTER_VALUE(OUT), type == 1 ? "wb" : "wt")))
		err(3, "Unable to open file `%s`\n", CHARACTER_VALUE(OUT));
	if(!(sp = rdfile(fin, type, swap, &c)))
		err(4, "Unable to read file `%s`\n", CHARACTER_VALUE(IN));
	if(!(wrfile(fout, sp, 1 - type, swap, header)))
		err(5, "Unable to write file `%s`\n", CHARACTER_VALUE(OUT));
	ret = 0;
error:
	spectra_free(sp);
	close(fin);
	close(fout);
	return mkInt(ret);
}
