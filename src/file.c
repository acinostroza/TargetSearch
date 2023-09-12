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

int checksig(unsigned char *s)
{
	int i = 0;
	for(i = 0; i < SIGLEN; i++)
		if(s[i] != signature[i])
			return 0;
	return 1;
}

/* from R source code: swap bytes for big endian machines */
void swapb(void *ptr, int size)
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
void swapp(void *ptr, int size, int len) {
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
void write_dat(FILE *fp, spectra_t *sp, int swap)
{
        int i, j, tmp, n, splen=sp->n_scans, pcount = 0;
	double ri, rt;

	/* write files signatures */
	fwrite(signature, SIGLEN, 1, fp);

        /* write len*/
	n = sp->n_scans;
	if(swap == 1)
		swapb(&n, sizeof(n));
        fwrite(&n, sizeof(int), 1, fp);

	/* write point count */
	for(i = 0; i < splen; i++)
		pcount += sp->sp[i].len;
	if(swap == 1)
		swapb(&pcount, sizeof(n));

        fwrite(&pcount, sizeof(int), 1, fp);

        /* write RI */
        for(i = 0; i < splen; i++) {
		ri = sp->ri[i];
		if(swap == 1)
			swapb(&ri, sizeof(ri));
                fwrite(&ri, sizeof(ri), 1, fp);
	}

        /* write RT */
        for(i = 0; i < splen; i++) {
		rt = sp->rt[i];
		if(swap == 1)
			swapb(&rt, sizeof(rt));
                fwrite(&rt, sizeof(rt), 1, fp);
	}

        /* write N */
        for(i = 0; i < splen; i++) {
		n = sp->sp[i].len;
		if(swap == 1)
			swapb(&n, sizeof(n));
                fwrite(&n, sizeof(n), 1, fp);
	}

        for(i = 0; i < splen; i++) {
		for(j = 0; j < sp->sp[i].len; j++) {
			tmp = sp->sp[i].mz[j];
			if(swap == 1)
				swapb(&tmp, sizeof(tmp));
			fwrite(&tmp, sizeof(tmp), 1, fp);
		}
		for(j = 0; j < sp->sp[i].len; j++) {
			tmp = sp->sp[i].in[j];
			if(swap == 1)
				swapb(&tmp, sizeof(tmp));
			fwrite(&tmp, sizeof(tmp), 1, fp);
		}
        }
}

/* macros needed to read a TXT file with peak data */
#define MAX(a, b)  ((a) > (b) ? (a) : (b))
#define CHKCOL(x, col) do \
	if((x) < 0) { \
		REprintf("Unable to find colum `%s'\n", col); \
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

/* R interface to write_dat functions */
void write_peaks_dat(char **fout, double *rt, double *ri, int *in, int *mass,
	int *nscans, int *swap)
{
	FILE *fp;
	spectra_t *sp = pktosp(rt, ri, in, mass, *nscans);
	if(sp == NULL)
		error("Error creating spectra struct.\n");

	fp = fopen(fout[0], "wb");
	if(fp == NULL)
		error("Error writing file %s\n", fout[0]);

	write_dat(fp, sp, *swap);
	spectra_free(sp);
	fclose(fp);
}

/* R interface to write_text functions */
void write_peaks_text(char **fout, double *rt, double *ri, int *in, int *mass,
	int *nscans, char **header)
{
	FILE *fp;
	spectra_t *sp = pktosp(rt, ri, in, mass, *nscans);
	if(sp == NULL)
		error("Error creating spectra struct.\n");

	fp = fopen(fout[0], "wb");
	if(fp == NULL)
		error("Error writing file %s\n", fout[0]);

	write_txt(fp, sp, header[0]);
	spectra_free(sp);
	fclose(fp);
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
 * @return void
 */
void write_txt(FILE *fp, spectra_t *sp, const char *header)
{
	fprintf(fp, "%s\n", header);
	for(int i = 0; i < sp->n_scans; i++) {
		spectrum_t *p = sp->sp + i;
		if(p->len == 0)
			continue;
		fprintf(fp, "%.15g\t", sp->rt[i]);
		for(int j = 0; j < p->len; j++)
			fprintf(fp, "%s%d:%d", j > 0 ? " " : "", p->mz[j], p->in[j]);
		fprintf(fp, "\t%.15g\n", sp->ri[i]);
	}
}

#define ERROR(...) do { snprintf(msg, sizeof(msg), __VA_ARGS__); goto clean; } while(0)
#define Fclose(fp) do { if(fp) fclose(fp); fp = NULL; } while(0)

/* Function to convert from TXT to DAT format*/
void text_to_dat(char **infile, char **outfile, int *swap, int *cols)
{
	FILE *fpin = NULL, *fpout = NULL;
	spectra_t *sp = NULL;
	// char *SPECTRUM_COL = cols[0], *RI_COL = cols[1], *RT_COL = cols[2];
	char msg[256] = {'\0'};
	int status = EXIT_FAILURE;

	if((fpin = fopen(infile[0], "r")) == NULL)
		ERROR("Error opening file %s\n", infile[0]);

	if((sp = read_txt(fpin, NULL, NULL, NULL, cols)) == NULL)
		ERROR("Error reading file %s\n", infile[0]);

	if((fpout = fopen(outfile[0], "wb")) == NULL)
		ERROR("Error opening file %s\n", outfile[0]);

	write_dat(fpout, sp, *swap);
	status = EXIT_SUCCESS;
clean:
	spectra_free(sp);
	Fclose(fpin);
	Fclose(fpout);
	if(status == EXIT_FAILURE)
		error("%s", msg);
}

/* Function to convert from DAT to TXT format*/
void dat_to_text(char **infile, char **outfile, int *swap, char **header)
{
	FILE *fpin = NULL, *fpout = NULL;
	char msg[256] = {'\0'};
	int status = EXIT_FAILURE;
	spectra_t *sp = NULL;

	if((fpin = fopen(infile[0], "rb")) == NULL)
		ERROR("Error opening file %s\n", infile[0]);

	if((sp = read_dat(fpin, *swap)) == NULL)
		ERROR("Error reading file %s\n", infile[0]);

	if((fpout = fopen(outfile[0], "w")) == NULL)
		ERROR("Error opening file %s\n", outfile[0]);

	write_txt(fpout, sp, header[0]);
	status = EXIT_SUCCESS;
clean:
	Fclose(fpin);
	Fclose(fpout);
	spectra_free(sp);
	if(status == EXIT_FAILURE)
		error("%s", msg);
}
