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
        {"FindPeaks", (DL_FUNC)&FindPeaks, 8},
        {"peak_finding", (DL_FUNC)&peak_finding, 7},
        {"ncdfToMatrix", (DL_FUNC)&ncdfToMatrix, 2},
        {"ppc", (DL_FUNC)&ppc, 4},
        {"cdffix", (DL_FUNC)&cdffix, 2},
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
