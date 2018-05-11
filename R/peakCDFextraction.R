`peakCDFextraction` <-
function(cdfFile, massRange = NULL)
{
	if(.is_ts_ncdf4(cdfFile))
		.peakCDFextraction4(cdfFile)
	else
		.peakCDFextraction3(cdfFile, massRange)
}

`.peakCDFextraction3` <-
function(cdfFile, massRange)
{
	ncData <- .open.ncdf(cdfFile)

	if(is.null(massRange))
		massRange <- range(ncData$mz)

	peaks <- .Call("ncdfToMatrix", ncData, massRange, PACKAGE="TargetSearch")
	colnames(peaks) <- as.character(massRange[1]:massRange[2])
	return( list(Time = ncData$rt, Peaks = peaks, massRange=massRange) )
}

`.peakCDFextraction4` <-
function(cdfFile, massRange)
{
	nc <- nc_open(cdfFile)
	on.exit(nc_close(nc))
	massRange <- ncvar_get(nc, 'mass_range')
	Time      <- ncvar_get(nc, 'retention_time')
	peaks     <- ncvar_get(nc, 'intensity')
	colnames(peaks) <- as.character(massRange[1]:massRange[2])
	time_corrected <-  ncatt_get(nc, 0, 'time_corrected')
	Index <- if(time_corrected$value==0) NULL else ncvar_get(nc, "retention_index")
	return(list(Time = Time, Peaks = peaks, massRange = massRange, Index=Index))
}

# vim: set ts=4 sw=4:
