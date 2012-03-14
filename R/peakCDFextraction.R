`peakCDFextraction` <-
function(cdfFile, massRange = c(85,500)) {

	nc     <- mzR:::netCDFOpen(cdfFile)
	ncData <- mzR:::netCDFRawData(nc)
	mzR:::netCDFClose(nc)
	point_count <- diff(c(ncData$scanindex, length(ncData$mz)))
	peaks <- .Call("peakExtraction", ncData$mz, ncData$intensity, point_count,
		ncData$scanindex, massRange, PACKAGE="TargetSearch" )
	return( list(Time = ncData$rt, Peaks = peaks) )
}

