# function to find peaks in chromatograms
# arguments:
#    - cdfFile : The file
#    - Window: Window to be used by method
#    - IntThreshold: Intensity threshold
#    - method: peak picking method. PPC / smoothing
#    - massRange: m/z range (deprecated, it has no effect)

NetCDFPeakFinding <- function(cdfFile, massRange = NULL, Window = 15, IntThreshold = 10,
	pp.method = "ppc", baseline = FALSE, baseline.opts = NULL) {

	method <- pmatch(pp.method,c("smoothing", "ppc", "gaussian"))

	if(is.na(method))
		stop("Invalid peak picking method.")

	method <- c("s", "p", "g")[method]

	refine <- min(c(round(Window / 3), 4))

	# first check that the file is not TS
	if(.is_ts_ncdf4(cdfFile)) {
		return(.PeakFinding(cdfFile, method, Window, IntThreshold, refine, baseline, baseline.opts))
	}

	ncData <- .open.ncdf(cdfFile)
	massRange <- range(ncData$mz)

	if(baseline)
		ncData <- baseline(ncData, baseline.opts)

	peaks <- .Call(c_peak_detection_main, method, ncData, Window, refine, IntThreshold, NULL)

	colnames(peaks) <- as.character(massRange[1]:massRange[2])
	return( list(Time = ncData$rt, Peaks = peaks, massRange=massRange) )
}

# peak finding for version 4
.PeakFinding <- function(cdfFile, method, Window, IntThreshold, refine, baseline,
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

	peaks <- .Call(c_peak_detection_main, method, NULL, Window, refine, IntThreshold, ncInt)

	colnames(peaks) <- as.character(mzRange[1]:mzRange[2])
	list(Time = Time, Peaks = peaks, massRange=mzRange)
}

# vim: set ts=4 sw=4:
