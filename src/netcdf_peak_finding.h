#ifndef _NETCDF_PEAK_FINDING_H
#define _NETCDF_PEAK_FINDING_H

SEXP peak_finding(SEXP MassValues, SEXP IntensityValues, SEXP PointCount,
	SEXP ScanIndex, SEXP Window, SEXP MassLimits, SEXP MinInt);

#endif
