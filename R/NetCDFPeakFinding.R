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

	ncData <- .open.ncdf(cdfFile)

	if(baseline)
		ncData <- baseline(ncData, baseline.opts)

	if(is.null(massRange))
		massRange <- range(ncData$mz)

	if(method == 1)
		peaks <- .Call("peak_finding", ncData$mz, ncData$intensity, ncData$point_count,
			ncData$scanindex, Window, massRange, IntThreshold, PACKAGE="TargetSearch")

	if(method == 2)
		peaks <- .Call("ppc", ncData, Window, massRange, IntThreshold, PACKAGE="TargetSearch")
	colnames(peaks) <- as.character(massRange[1]:massRange[2])
	return( list(Time = ncData$rt, Peaks = peaks, massRange=massRange) )
}

# vim: set ts=4 sw=4:
