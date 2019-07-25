#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* .Call interface */
/* FindPeaks */
#include "find_peaks.h"
/* cdffix and ncdf_to_matrix */
#include "ncdf.h"
/* from dectection.c */
SEXP peak_detection_main(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* from baseline.c */
SEXP baseline(SEXP, SEXP, SEXP, SEXP, SEXP);

/* .C interface */
/* hpf, windowing for baseline correction */
#include "hpf.h"
/* write_peaks_dat, write_peaks_text, text_to_dat, dat_to_text */
#include "file.h"

/* Automate using sed or something. */
#if _MSC_VER >= 1000
__declspec(dllexport)
#endif

static const R_CallMethodDef R_CallDef[] = {
        {"find_peaks", (DL_FUNC)&find_peaks, 8},
        {"ncdf_to_matrix", (DL_FUNC)&ncdf_to_matrix, 2},
        {"nominal", (DL_FUNC)&nominal, 1},
        {"peak_detection_main", (DL_FUNC)&peak_detection_main, 6},
        {"baseline", (DL_FUNC)&baseline, 5},
        {NULL, NULL, 0},
};

static R_NativePrimitiveArgType hpf_t[] = {REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType windowing_t[] = {
                INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType write_peaks_dat_t[] = {
                STRSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType write_peaks_text_t[] = {
                STRSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, STRSXP};
static R_NativePrimitiveArgType text_to_dat_t[] = {STRSXP, STRSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType dat_to_text_t[] = {STRSXP, STRSXP, INTSXP, STRSXP};

static const R_CMethodDef cMethods[] = {
        {"hpf", (DL_FUNC)&hpf, 4, hpf_t},
        {"windowing", (DL_FUNC)&windowing, 5, windowing_t},
        {"write_peaks_dat", (DL_FUNC)&write_peaks_dat, 7, write_peaks_dat_t},
        {"write_peaks_text", (DL_FUNC)&write_peaks_text, 7, write_peaks_text_t},
        {"text_to_dat", (DL_FUNC)&text_to_dat, 4, text_to_dat_t},
        {"dat_to_text", (DL_FUNC)&dat_to_text, 4, dat_to_text_t},
        {NULL, NULL, 0, NULL}
};

void R_init_TargetSearch(DllInfo *info)
{
        R_registerRoutines(info, cMethods, R_CallDef, NULL, NULL);
        R_useDynamicSymbols(info, FALSE);
}
