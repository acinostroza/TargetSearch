#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/RConverters.h>
#include <R_ext/Rdynload.h>

SEXP FindPeaks(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP peak_finding(SEXP, SEXP, SEXP,	SEXP, SEXP, SEXP, SEXP);
SEXP peakExtraction(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP ppc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* Automate using sed or something. */
#if _MSC_VER >= 1000
__declspec(dllexport)
#endif
	
static const R_CallMethodDef R_CallDef[] = {
        {"FindPeaks", (DL_FUNC)&FindPeaks, 5},
        {"peak_finding", (DL_FUNC)&peak_finding, 7},
        {"peakExtraction", (DL_FUNC)&peakExtraction, 5},
        {"ppc", (DL_FUNC)&ppc, 7},
        {NULL, NULL, 0},
};

void R_init_TargetSearch(DllInfo *info)
{
  R_registerRoutines(info,NULL,R_CallDef,NULL,NULL);
}
