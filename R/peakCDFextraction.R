`peakCDFextraction` <-
function(cdfFile, massRange = NULL)
{
	ncData <- .open.ncdf(cdfFile)

	if(is.null(massRange))
		massRange <- range(ncData$mz)

	peaks <- .Call("ncdfToMatrix", ncData, massRange, PACKAGE="TargetSearch")
	colnames(peaks) <- as.character(massRange[1]:massRange[2])
	return( list(Time = ncData$rt, Peaks = peaks, massRange=massRange) )
}

# vim: set ts=4 sw=4:
