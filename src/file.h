#ifndef _FILE_H
#define _FILE_H

#define STRLEN 256

#include "spectrum.h"

/* helper struct for parsing column options from R */
struct column_s {
	const char *sp_col;
	const char *ri_col;
	const char *rt_col;
	int * icols;
};

/* Function prototypes */
spectra_t * read_dat(FILE *, int);
int write_dat(FILE *, spectra_t *, int);

spectra_t * read_txt(FILE *fp, const char *, const char *, const char *, const int *);
int write_txt(FILE *, spectra_t *, const char *);

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
SEXP write_peaks(SEXP Output, SEXP RT, SEXP RI, SEXP IN, SEXP Mass,
		SEXP Swap, SEXP Header);

/**
 * C interface to write RI files from and to text to bin
 *
 * @param IN (string) path to the input file to read.
 * @param OUT (string) path to the output file to write.
 * @param Type (integer) whether to convert bin2text (0) or text2bin (1).
 * @param Columns the column names or indices of the text file.
 * @param Header the header of the tab-delimited file to write (string).
 *
 * The columns should be in the order 'spectrum', 'ret_index', 'ret_time'
 * while the header should be 'ret_time', 'spectrum', 'ret_index'. This
 * is due to historical reasons.
 *
 * @return An integer. A non-zero value means an error occurred, otherwise
 *        returns zero.
 */
SEXP convert_ri_file(SEXP IN, SEXP OUT, SEXP Type, SEXP Columns, SEXP Header);

#endif
