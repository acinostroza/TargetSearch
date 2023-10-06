#ifndef SPECTRUM_H
#define SPECTRUM_H

struct _spectrum {
	int *mz;   /* mass or m/z */
	int *in;   /* intensity */
	int len;   /* number of points */
	int alloc; /* allocated memory */
};

typedef struct _spectrum spectrum_t;

struct _spectra {
	int n_scans;  /* number of scans */
	int alloc;    /* size of allocated memory */
	double *ri;
	double *rt;
	struct _spectrum *sp;
};

typedef struct _spectra spectra_t;

int spectrum_init(spectrum_t *z, int alloc);
int spectrum_add(spectrum_t *z, int mz, int in);
spectra_t * spectra_init(int alloc);
int spectra_add(spectra_t * spectra, double rt, double ri, spectrum_t *z);

void _spectra_free(spectra_t * spectra);
#define spectra_free(p) do { _spectra_free(p); (p) = NULL; } while(0)

/* error codes for spectrum_parse */
#define SP_ERR_INVALID  -1
#define SP_ERR_STRTOL   -2
#define SP_ERR_OVERFLOW -3
#define SP_ERR_ALLOC    -4

int spectrum_parse(const char *str, spectrum_t * z);
int spectrum_scan(const char *str);

#endif
