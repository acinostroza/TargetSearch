#ifndef _FILE_H
#define _FILE_H

#define STRLEN 256

#include "spectrum.h"

/* Function prototypes */
spectra_t * read_dat(FILE *, int);
void write_dat(FILE *, spectra_t *, int);
void write_peaks_dat(char **, double *, double *, int *, int *, int *, int *);

spectra_t * read_txt(FILE *fp, const char *, const char *, const char *, const int *);
void write_txt(FILE *, spectra_t *, const char *);
void write_peaks_text(char **, double *, double *, int *, int *, int *, char **);

spectra_t * pktosp(double *, double *, int *, int *, int);

void text_to_dat(char **, char **, int *, int *);
void dat_to_text(char **, char **, int *, char **);

#endif
