`peakCDFextraction` <-
function(cdfFile, massRange = NULL)
{
	nc     <- mzR:::netCDFOpen(cdfFile)
	ncData <- mzR:::netCDFRawData(nc)
	mzR:::netCDFClose(nc)
	ncData$point_count <- diff(c(ncData$scanindex, length(ncData$mz)))
	ncData <- .check.mz.precision(ncData)
	if(is.null(ncData)) {
		stop(paste("Error processing file '", cdfFile, "'. It seems to contains",
			"high mass accuracy data.", sep=" "))
	}

	if(is.null(massRange))
		massRange <- range(ncData$mz)
	peaks <- .Call("peakExtraction", ncData$mz, ncData$intensity,
		ncData$point_count, ncData$scanindex, massRange,
		PACKAGE="TargetSearch")
	colnames(peaks) <- as.character(massRange[1]:massRange[2])
	return( list(Time = ncData$rt, Peaks = peaks, massRange=massRange) )
}

# vim: set ts=4 sw=4:
