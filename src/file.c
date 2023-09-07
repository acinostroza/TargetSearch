/* Functions to perform file manipulation */
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <ctype.h>
#include "file.h"
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
 * returns: spectra_t struct */
spectra_t * read_dat(FILE *fp, int swap)
{
	spectra_t * spectra;
	int i, splen=0, sptotal=0, *pcount = NULL;
	unsigned char s[SIGLEN];

	spectra = (spectra_t *) R_alloc(1, sizeof(spectra_t));

	if(!fread(s, SIGLEN, 1, fp))
		return NULL;
	if(checksig(s) == 0) {
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

	spectra->n_scans = splen;
	spectra->p_count = sptotal;

	spectra->ri = (double *) R_alloc(splen, sizeof(double));
	spectra->rt = (double *) R_alloc(splen, sizeof(double));
	pcount  = (int *) R_alloc(splen, sizeof(int));
	spectra->sp = (spectrum_t *) R_alloc(splen, sizeof(spectrum_t));

	if(!fread(spectra->ri, splen * sizeof(double), 1, fp))
		return NULL;
	if(!fread(spectra->rt, splen * sizeof(double), 1, fp))
		return NULL;
	if(!fread(pcount, splen *  sizeof(int), 1, fp))
		return NULL;

	if(swap == 1) {
		swapp(pcount, sizeof(int), splen);
		swapp(spectra->ri, sizeof(double), splen);
		swapp(spectra->rt, sizeof(double), splen);
	}

	spectra->sp[0].mz   = (int *) R_alloc(sptotal, sizeof(int));
	spectra->sp[0].in   = (int *) R_alloc(sptotal, sizeof(int));

	for(i = 0; i < splen; i++) {
		spectrum_t * z = spectra->sp + i;
		z->len = pcount[i];
		if(i > 0) {
			spectra->sp[i].mz   = spectra->sp[i-1].mz   + pcount[i-1];
			spectra->sp[i].in   = spectra->sp[i-1].in   + pcount[i-1];
		}
		if(!fread(spectra->sp[i].mz, pcount[i] * sizeof(int), 1, fp))
			return NULL;
		if(!fread(spectra->sp[i].in, pcount[i] * sizeof(int), 1, fp))
			return NULL;
		if(swap == 1) {
			swapp(spectra->sp[i].mz, sizeof(int), pcount[i]);
			swapp(spectra->sp[i].in, sizeof(int), pcount[i]);
		}
	}
	return spectra;
}

/* function to write a dat file */
void write_dat(FILE *fp, spectra_t *sp, int swap)
{
        int i, j, tmp, n, splen=sp->n_scans;
	double ri, rt;

	/* write files signatures */
	fwrite(signature, SIGLEN, 1, fp);

        /* write len*/
	n = sp->n_scans;
	if(swap == 1)
		swapb(&n, sizeof(n));
        fwrite(&n, sizeof(int), 1, fp);

	/* write point count */
	n = sp->p_count;
	if(swap == 1)
		swapb(&n, sizeof(n));
        fwrite(&n, sizeof(int), 1, fp);

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

/* function to read a spectra file in TAB-DELIMITED format */
spectra_t * read_txt(FILE *fp, int SPECTRUM_COL, int RI_COL, int RT_COL)
{
	spectra_t *spectra;
	int total = 0, j;
	int header = 1;
	char *line = NULL;
	int  len = 0;
	int  err = 0;
	int  next = EMPTY_CHAR;

	char *ri_str = NULL, *sp_str = NULL, *rt_str = NULL;
	int  ri_i, sp_i, rt_i, tabs, n;
	int  ri_len = 0, sp_len = 0, rt_len = 0;

	spectra = (spectra_t *) R_alloc(1, sizeof(spectra_t));

	while (get_line(&line, &len, &next, fp) > 0) {
		total++;
	}

	total--; /* header isn't counted */
	spectra->n_scans = total;
	spectra->p_count = 0;

	spectra->ri = (double *) R_alloc(total, sizeof(double));
	spectra->rt = (double *) R_alloc(total, sizeof(double));
	spectra->sp = (spectrum_t *) R_alloc(total, sizeof(spectrum_t));

	fseek(fp, 0, SEEK_SET);
	j = 0;

	while (get_line(&line, &len, &next, fp) > 0) {
		size_t i = 0;
		if (header) {
			header = 0;
			continue;
		}
		tabs = ri_i = rt_i = sp_i = n = 0;

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

		if(n == 0 || ri_i == 0 || rt_i == 0 || sp_i == 0) {
			REprintf("Error reading spectra. Invalid spectrum format:\n");
			REprintf("--> Line %d: '%s'\n", j+1, line);
			err = 1;
			goto end;
		}
		spectra->p_count += n;

		spectra->sp[j].len = n;
		spectra->ri[j] = atof(ri_str);
		spectra->rt[j] = atof(rt_str);

		spectra->sp[j].mz = (int *) R_alloc(n , sizeof(int));
		spectra->sp[j].in = (int *) R_alloc(n , sizeof(int));

		if(read_spectrum(sp_str, spectra->sp[j].mz, spectra->sp[j].in) == 0) {
			REprintf("Error reading spectra. Invalid spectrum format:\n");
			REprintf("--> Line %d: '%s'\n", j+1, line);
			err = 1;
			goto end;
		}
		j++;
	}

end:
	if(line)
		R_chk_free(line);
	if(ri_str)
		R_chk_free(ri_str);
	if(rt_str)
		R_chk_free(rt_str);
	if(sp_str)
		R_chk_free(sp_str);
	if(err)
		return NULL;
	return spectra;
}

/* function to parse a spectrum string *
 * the format is a list of mass:intensity pairs separated by spaces *
 * spectrum: the spectrum string.
 * mass, in; the returned masses and intensities vectors.
*/
int read_spectrum(char *spectrum, int *mass, int *in)
{
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

/* transform a peak matrix to spectra_t struct */
spectra_t pktosp(double *rt, double *ri, int *in, int *mass, int nscans)
{
	spectra_t sp;
	int i, j, k=0, l=0, splen=0, sptotal=0, n[nscans];
	int nmz=mass[1]-mass[0]+1;

	sp.n_scans = -1;
	sp.sp = NULL;
	sp.rt = NULL;
	sp.ri = NULL;
	sp.p_count = 0;

	if(nmz < 0)
		return sp;

	/* scan data */
	for(i = 0; i < nscans; i++) {
		n[i] = 0;
		for(j = 0; j < nmz; j++) {
			if(in[i + j * nscans] > 0)
				n[i]++;
		}
		if(n[i] > 0) {
			splen++;
			sptotal += n[i];
		}
	}
	sp.n_scans = splen;
	sp.p_count = sptotal;

	sp.ri         = (double *) R_alloc(splen, sizeof(double));
	sp.rt         = (double *) R_alloc(splen, sizeof(double));
	sp.sp         = (spectrum_t *) R_alloc(splen, sizeof(spectrum_t));
	sp.sp[0].mz   = (int *) R_alloc(sptotal, sizeof(int));
	sp.sp[0].in   = (int *) R_alloc(sptotal, sizeof(int));

	/* construct the structure */
	for(i = 0; i < nscans; i++) {
		if(n[i] == 0)
			continue;
		sp.sp[k].len = n[i];
		sp.rt[k] = rt[i];
		sp.ri[k] = ri[i];
		if(k > 0) {
			sp.sp[k].mz = sp.sp[k-1].mz + sp.sp[k - 1].len;
			sp.sp[k].in = sp.sp[k-1].in + sp.sp[k - 1].len;
		}
		l = 0;
		for(j = 0; j < nmz; j++) {
			if(in[i + j * nscans] > 0) {
				sp.sp[k].mz[l] = mass[0] + j;
				sp.sp[k].in[l++] = in[i + j * nscans];
			}
		}
		k++;
	}
	return sp;
}

/* R interface to write_dat function */
void write_peaks_dat(char **fout, double *rt, double *ri, int *in, int *mass,
	int *nscans, int *swap)
{
	FILE *fp;
	spectra_t sp;

	sp = pktosp(rt, ri, in, mass, *nscans);
	if(sp.n_scans == -1)
		error("Error creacting spectra struct.\n");

	fp = fopen(fout[0], "wb");
	if(fp == NULL)
		error("Error writing file %s\n", fout[0]);

	write_dat(fp, &sp, *swap);
	fclose(fp);
}

/* R interface to write_txt function */
void write_peaks_text(char **fout, double *rt, double *ri, int *in, int *mass,
	int *nscans, char **header)
{
	FILE *fp;
	spectra_t sp;

	sp = pktosp(rt, ri, in, mass, *nscans);
	if(sp.n_scans == -1)
		error("Error creacting spectra struct\n");

	fp = fopen(fout[0], "w");
	if(fp == NULL)
		error("Error writing file %s\n", fout[0]);

	write_txt(fp, &sp, header[0]);
	fclose(fp);
}

void write_txt(FILE *fp, spectra_t *sp, char *header)
{
	int i, j;
	char c;
	fprintf(fp, "%s\n", header);
	for(i = 0; i < sp->n_scans; i++) {
		fprintf(fp, "%.15g\t", sp->rt[i]);
		for(j = 0; j < sp->sp[i].len; j++) {
			fprintf(fp, "%d:%d", sp->sp[i].mz[j], sp->sp[i].in[j]);
			c = (j == sp->sp[i].len - 1) ? '\t' : ' ';
			fputc(c, fp);
		}
		fprintf(fp, "%.15g\n", sp->ri[i]);
	}
}

/* Function to convert from TXT to DAT format*/
void text_to_dat(char **infile, char **outfile, int *swap, int *cols)
{
	FILE *fpin, *fpout;
	spectra_t *sp;

	fpin = fopen(infile[0], "r");
	if(fpin == NULL)
		error("Error opening file %s\n", infile[0]);
	fpout = fopen(outfile[0], "wb");
	if(fpout == NULL)
		error("Error opening file %s\n", outfile[0]);
	sp = read_txt(fpin, cols[0], cols[1], cols[2]);
	if(!sp)
		error("Error reading file %s\n", infile[0]);
	write_dat(fpout, sp, *swap);
	fclose(fpin);
	fclose(fpout);
}

/* Function to convert from DAT to TXT format*/
void dat_to_text(char **infile, char **outfile, int *swap, char **header)
{
	FILE *fpin, *fpout;
	spectra_t *sp;

	fpin = fopen(infile[0], "rb");
	if(fpin == NULL)
		error("Error opening file %s\n", infile[0]);
	fpout = fopen(outfile[0], "w");
	if(fpout == NULL)
		error("Error opening file %s\n", outfile[0]);
	sp = read_dat(fpin, *swap);
	if(!sp)
		error("Error reading file %s\n", infile[0]);
	write_txt(fpout, sp, header[0]);
	fclose(fpin);
	fclose(fpout);
}
