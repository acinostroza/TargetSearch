# function to find peaks in chromatograms
# arguments:
#    - cdfFile : The file
#    - massRange: m/z range.
#    - Window: Window to be used by method
#    - IntThreshold: Intensity threshold
#    - method: peak picking method. PPC / smoothing

NetCDFPeakFinding <- function(cdfFile, massRange = NULL, Window = 15, IntThreshold = 10,
	pp.method = "ppc", baseline = FALSE, baseline.opts = NULL) {

	method <- pmatch(pp.method,c("smoothing", "ppc"))

	if(is.na(method))
		stop("Invalid peak picking method.")

	# first check that the file is not TS
	if(.is_ts_ncdf4(cdfFile)) {
		return(.PeakFinding(cdfFile, Window, IntThreshold, method, baseline, baseline.opts))
	}

	ncData <- .open.ncdf(cdfFile)
	if(baseline)
		ncData <- baseline(ncData, baseline.opts)

	if(is.null(massRange))
		massRange <- range(ncData$mz)

	if(method == 1)
		peaks <- .Call("peak_finding", ncData$mz, ncData$intensity, ncData$point_count,
			ncData$scanindex, Window, massRange, IntThreshold, PACKAGE="TargetSearch")

	if(method == 2)
		peaks <- .Call("ppc", ncData, Window, massRange, IntThreshold, NULL, PACKAGE="TargetSearch")
	colnames(peaks) <- as.character(massRange[1]:massRange[2])
	return( list(Time = ncData$rt, Peaks = peaks, massRange=massRange) )
}

# peak finding for version 4
.PeakFinding <- function(cdfFile, Window, IntThreshold, method, baseline,
	baseline.opts = NULL)
{
	nc      <- nc_open(cdfFile)
	ncInt   <- ncvar_get(nc, 'intensity')
	Time    <- ncvar_get(nc, 'retention_time')
	mzRange <- ncvar_get(nc, 'mass_range')
	nc_close(nc)

	if(baseline) {
		ncInt <- do.call(baselineCorrection, append(list(int = t(ncInt)), baseline.opts))
		ncInt <- t(ncInt)
	}

	if(method == 1)
		stop('Error: method not implemented for netCDF-4. use a netCDF-3 file')

	if(method == 2)
		peaks <- .Call("ppc", NULL, Window, NULL, IntThreshold, ncInt, PACKAGE="TargetSearch")

	colnames(peaks) <- as.character(mzRange[1]:mzRange[2])
	list(Time = Time, Peaks = peaks, massRange=mzRange)
}

# vim: set ts=4 sw=4:
