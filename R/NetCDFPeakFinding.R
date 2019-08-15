# function to find peaks in chromatograms
# arguments:
#    - ncdf : The file (version 3 or 4) or a alist
#    - Window: Window to be used by method
#    - IntThreshold: Intensity threshold
#    - method: peak picking method. PPC / smoothing
#    - massRange: m/z range (deprecated, it has no effect)

NetCDFPeakFinding <- function(ncdf, massRange = NULL, Window = 15, IntThreshold = 10,
	pp.method = "ppc", baseline = FALSE, ...) {

	method <- match.arg(pp.method,c("smoothing", "ppc", "gaussian"))
	method <- substring(method, 1, 1)

	# methods smoothing and gaussian use W points, while ppc uses 2*W + 1
	# make the number of points homogeneous,
	Window <- switch(method, s=2*Window+1, p=Window, g=2*Window+1)

	refine <- min(c(round(Window / 3), 4))

	ncData <- .peakExtractWrap(ncdf)

	if(baseline)
		ncData <- baseline(ncData, ...)

	peaks <- .Call(c_peak_detection_main, method, NULL, Window, refine, IntThreshold, ncData$Peaks)

	massRange <- ncData$massRange
	colnames(peaks) <- as.character(massRange[1]:massRange[2])
	return( list(Time = ncData$Time, Peaks = peaks, massRange=massRange) )
}

# vim: set ts=4 sw=4:
