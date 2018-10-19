#ifndef _FILE_H
#define _FILE_H

#define is_open(x,f,p)  if (x == NULL) { \
			UNPROTECT(p); \
			error("Unable to open file %s.\n", f); \
			return(R_NilValue); \
			}

/* memory allocation macro */
#define str_alloc(x,x_len,len) if(len > x_len) { \
                        	if(x) \
                                	R_chk_free(x); \
	                        x = (char *) R_chk_calloc(len, sizeof(char)); \
	                        x_len = len; \
                		} 

#define STRLEN 256

struct peaks {
	int  *mass; /* mass */
	int    *in; /* intensity */
};

typedef struct {
	int n_scans;  /* number of scans */
	int p_count;  /* point count */
	double *ri;
	double *rt;
	int *n;
	struct peaks *pk;
} SPECTRA;

/* Function prototypes */
int checksig(unsigned char *);

void swapb(void *, int);
void swapp(void *, int, int);

SPECTRA * read_dat(FILE *, int);
void write_dat(FILE *, SPECTRA *, int);
void write_peaks_dat(char **, double *, double *, int *, int *, int *, int *);

SPECTRA * read_txt(FILE *, int, int, int);
int read_spectrum(char *, int *, int *);
void write_txt(FILE *, SPECTRA *, char *);
void write_peaks_text(char **, double *, double *, int *, int *, int *, char **);

SPECTRA pktosp(double *, double *, int *, int *, int);

void text_to_dat(char **, char **, int *, int *);
void dat_to_text(char **, char **, int *, char **);

#endif
