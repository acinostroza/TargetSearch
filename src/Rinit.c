#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

SEXP FindPeaks(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP peak_finding(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP peakExtraction(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP ppc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP cdffix(SEXP, SEXP, SEXP, SEXP);

void hpf(double *, double *, int *, double *);
void windowing(int *, int *, int *, int *, int *);
void writePeaksDAT(char **, double *, double *, int *, int *, int *, int *);
void writePeaksTXT(char **, double *, double *, int *, int *, int *, char **);
void txt2dat(char **, char **, int *, int *);
void dat2txt(char **, char **, int *, char **);

/* Automate using sed or something. */
#if _MSC_VER >= 1000
__declspec(dllexport)
#endif

static const R_CallMethodDef R_CallDef[] = {
        {"FindPeaks", (DL_FUNC)&FindPeaks, 6},
        {"peak_finding", (DL_FUNC)&peak_finding, 7},
        {"peakExtraction", (DL_FUNC)&peakExtraction, 5},
        {"ppc", (DL_FUNC)&ppc, 7},
        {"cdffix", (DL_FUNC)&cdffix, 4},
        {NULL, NULL, 0},
};

static const R_CMethodDef cMethods[] = {
        {"hpf", (DL_FUNC)&hpf, 4},
        {"windowing", (DL_FUNC)&windowing, 5},
        {"writePeaksDAT", (DL_FUNC)&writePeaksDAT, 7},
        {"writePeaksTXT", (DL_FUNC)&writePeaksTXT, 7},
        {"txt2dat", (DL_FUNC)&txt2dat, 4},
        {"dat2txt", (DL_FUNC)&dat2txt, 4},
        {NULL, NULL, 0}
};

void R_init_TargetSearch(DllInfo *info)
{
  R_registerRoutines(info,cMethods,R_CallDef,NULL,NULL);
}
