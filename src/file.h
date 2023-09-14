#ifndef _FILE_H
#define _FILE_H

#define STRLEN 256

#include "spectrum.h"

/* Function prototypes */
spectra_t * read_dat(FILE *, int);
void write_dat(FILE *, spectra_t *, int);

spectra_t * read_txt(FILE *fp, const char *, const char *, const char *, const int *);
void write_txt(FILE *, spectra_t *, const char *);

spectra_t * pktosp(double *, double *, int *, int *, int);

void text_to_dat(char **, char **, int *, int *);
void dat_to_text(char **, char **, int *, char **);

/*
 * C interface to write RI files (peak list) from R
 *
 * @param Output path to the file name to write (string)
 * @param RT vector of retention time values (double)
 * @param RI vector of retention indices (double)
 * @param IN matrix of peak intensities (integer). A peak correspond with
 *        non-zero values.
 * @param Mass integer vector representing the m/z range (min, max).
 * @param Swap an integer. If -1, then it writes a text file; if 0, then
 *        write a binary file using little endian (no swapping); if 1, then
 *        perform swappin and store in little endian regardless.
 * @param Header string. The file header. Only used for writing text files,
 *        ignored otherwise.
 * @return TRUE on success, FALSE otherwise.
 */
SEXP write_peaks(SEXP Output, SEXP RT, SEXP RI,	SEXP IN, SEXP Mass,
		SEXP Swap, SEXP Header);

#endif
