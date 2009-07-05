# function to find peaks in chromatograms
# arguments: 
#    - cdfFile : The file
#    - massRange: m/z range.
#    - Window: Window to be used by method
#    - IntThreshold: Intensity threshold
#    - method: peak picking method. PPC / smoothing

NetCDFPeakFinding <- function(cdfFile, massRange =c(85,500), Window = 5, IntThreshold = 10,
	pp.method = "smoothing", baseline = FALSE, baseline.opts = NULL) {

	method <- pmatch(pp.method,c("smoothing", "ppc"))
	 
	if(is.na(method))
	  stop("Invalid peak picking method.")

	require(xcms)
	nc     <- xcms:::netCDFOpen(cdfFile)
	ncData <- xcms:::netCDFRawData(nc)
	xcms:::netCDFClose(nc)

	ncData$point_count <- diff(c(ncData$scanindex, length(ncData$mz)))
	
	if(baseline)
		ncData <- baseline(ncData, baseline.opts)

	if(method == 1)
		peaks <- .Call("peak_finding", ncData$mz, ncData$intensity, ncData$point_count,
			ncData$scanindex, Window, massRange, IntThreshold, PACKAGE="TargetSearch")

	if(method == 2)
		peaks <- .Call("ppc",  ncData$mz, ncData$intensity, ncData$point_count,
			ncData$scanindex, Window, massRange, IntThreshold, PACKAGE="TargetSearch")

	return( list(Time = ncData$rt, Peaks = peaks) )
}

