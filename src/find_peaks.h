
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

typedef struct {
	int      n; /* spectrum length */
	double  ri; /* retention time index */
	int  *mass; /* mass */
	int    *in; /* intensity */
} SPECTRA;

/* Function prototypes */

SPECTRA *read_file(FILE *fp, int *total, int, int);

int read_spectrum(char *spectrum, int *mass, int *in, int n);
void find_peak (double ,int ,double ,SPECTRA *,int ,double *, int *);

