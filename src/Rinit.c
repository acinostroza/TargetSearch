#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* FindPeaks */
#include "find_peaks.h"
/* peak_finding */
#include "netcdf_peak_finding.h"
/* ppc */
#include "ppc.h"
/* cdffix */
#include "ncdf.h"

/* hpf, windowing */
#include "hpf.h"
/* writePeaksDAT, writePeaksTXT, txt2dat, dat2txt */
#include "file.h"

/* Automate using sed or something. */
#if _MSC_VER >= 1000
__declspec(dllexport)
#endif

static const R_CallMethodDef R_CallDef[] = {
        {"find_peaks", (DL_FUNC)&find_peaks, 8},
        {"peak_finding", (DL_FUNC)&peak_finding, 7},
        {"ncdf_to_matrix", (DL_FUNC)&ncdf_to_matrix, 2},
        {"ppc", (DL_FUNC)&ppc, 5},
        {"cdffix", (DL_FUNC)&cdffix, 2},
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
