`peakCDFextraction` <-
function(cdfFile, massRange = NULL)
{
	ncData <- .open.ncdf(cdfFile)

	if(any(ncData$scanindex < 0)) {
		message('Warning:')
		message('  The following file seems to be corrupted. TargetSearch will attempt to process it anyway...')
		message(paste('  ->', cdfFile))
		# removing negative values
		tmp <- ncData$scanindex >= 0
		ncData$scanindex   <- ncData$scanindex[tmp]
		ncData$rt          <- ncData$rt[tmp]
		ncData$point_count <- ncData$point_count[tmp]
	}

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
