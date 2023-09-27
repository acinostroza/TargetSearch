#include <Rinternals.h>
#include <Rdefines.h>

/* .Call interface */
/* FindPeaks */
#include "find.h"
/* cdffix and ncdf_to_matrix */
#include "ncdf.h"
/* from dectection.c */
SEXP peak_detection_main(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* from baseline.c */
SEXP baseline(SEXP, SEXP, SEXP, SEXP, SEXP);

/* .C interface */
/* hpf, windowing for baseline correction */
#include "hpf.h"
/* write_peaks, text_to_dat, dat_to_text */
#include "file.h"

static const R_CallMethodDef R_CallDef[] = {
        {"find_peaks", (DL_FUNC)&find_peaks, 8},
        {"ncdf_to_matrix", (DL_FUNC)&ncdf_to_matrix, 2},
        {"nominal", (DL_FUNC)&nominal, 1},
        {"peak_detection_main", (DL_FUNC)&peak_detection_main, 6},
        {"baseline", (DL_FUNC)&baseline, 5},
        {"write_peaks", (DL_FUNC)&write_peaks, 7},
        {"convert_ri_file", (DL_FUNC)&convert_ri_file, 5},
        {NULL, NULL, 0},
};

static R_NativePrimitiveArgType hpf_t[] = {REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType windowing_t[] = {
                INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};

static const R_CMethodDef cMethods[] = {
        {"hpf", (DL_FUNC)&hpf, 4, hpf_t},
        {"windowing", (DL_FUNC)&windowing, 5, windowing_t},
        {NULL, NULL, 0, NULL}
};

void R_init_TargetSearch(DllInfo *info)
{
        R_registerRoutines(info, cMethods, R_CallDef, NULL, NULL);
        R_useDynamicSymbols(info, FALSE);
}

/* vim: set et: */
