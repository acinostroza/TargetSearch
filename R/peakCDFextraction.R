# peak CDF extraction converts a cdf into a matrix.
# The massRange is deprecated, so we ignore it
`peakCDFextraction` <-
function(cdfFile, massRange)
{
	if(.is_ts_ncdf4(cdfFile))
		.peakCDFextraction4(cdfFile)
	else
		.peakCDFextraction3(cdfFile)
}

`.peakCDFextraction3` <-
function(cdfFile)
{
	ncData <- .open.ncdf(cdfFile)
	massRange <- range(ncData$mz)

	peaks <- .Call(c_ncdf_to_matrix, ncData, massRange, PACKAGE="TargetSearch")
	colnames(peaks) <- as.character(massRange[1]:massRange[2])
	return( list(Time = ncData$rt, Peaks = peaks, massRange=massRange, Index=NULL,
				 baselineCorrected=FALSE) )
}

`.peakCDFextraction4` <-
function(cdfFile)
{
	nc <- nc_open(cdfFile)
	on.exit(nc_close(nc))
	massRange <- ncvar_get(nc, 'mass_range')
	Time      <- ncvar_get(nc, 'retention_time')
	peaks     <- ncvar_get(nc, 'intensity')
	colnames(peaks) <- as.character(massRange[1]:massRange[2])
	time_corrected <-  ncatt_get(nc, 0, 'time_corrected')
	Index <- if(time_corrected$value==0) NULL else ncvar_get(nc, "retention_index")

	basecor <- ncatt_get(nc, 0, 'baseline_corrected')
	basecor <- if(basecor$hasatt) basecor$value else 0

	return(list(Time = Time, Peaks = peaks, massRange = massRange, Index=Index,
				baselineCorrected=as.logical(basecor)))
}

#' Validates a NCDF strict (better with S4 class)
.validNCDF <- function(x) {
	if(!is.list(x))
		stop("Object must be a list")
	if(is.null(x$Time))
		stop("Missing `Time` attribute")
	if(is.null(x$Peaks))
		stop("Missing `Peaks` attribute")
	if(is.null(x$massRange))
		stop("Missing `massRange` attribute")
	nScan <- length(x$Time)
	nInd <- if(!is.null(x$Index)) length(x$Index) else nScan
	if(x$massRange[2] - x$massRange[1] + 1 != ncol(x$Peaks))
		stop("Invalid `massRange` or `Peaks` attribute. Not equal dims")
	if(nScan != nInd)
		stop("Invalid length of scans and indices")
	if(nrow(x$Peaks) != nScan)
		stop("Invalid length of scans and peak data")
	return(TRUE)
}

#'
#' Open CDF file and extract peaks or check CDF data.
#'
#' Detects whether the input is a string, in which case it assumes
#' is a CDF file and extracts it. If is a list, then it checks that
#' it is a valid CDF structure. In any other case, throws an error.
#'
#' @param ncdf a file name, ncdf list or matrix
#' @return a ncdf list struct.
`.peakExtractWrap` <- function(ncdf)
{
	if(is.character(ncdf))
		return(peakCDFextraction(ncdf[1]))
	if(is.list(ncdf)) {
		.validNCDF(ncdf)
		return(ncdf)
	}
	stop("Invalid parameter `ncdf`. Expecting a list or a file name")
}

# vim: set ts=4 sw=4:
