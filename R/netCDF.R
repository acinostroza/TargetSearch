# auxiliary functions to NetCDFPeakFinding and peakCDFextraction

# test whether the mass values are integers of double.
# if double, then checks whether they correspond with high mass
# accuracy data or not. If not, the mz values are rounded
# to the closes integer value.

.check.mz.precision <- function(x, force=FALSE) {
	mz <- sample(x$mz, min(c(5000, length(x$mz))))
	z <- sum(abs(mz - round(mz)))
	w <- all.equal(z, 0L)

	if(is.logical(w)) {
		if(w) {
			return(x)
		}
	}

	if(force) {
		x$mz <- as.integer(round(x$mz))
		return(.cdffix(x, max_assigned=100000L))
	}

	k <- which.max(x$point_count)
	mz <- x$mz[ x$scanindex[k] + c(1:x$point_count[k])]
	z  <- sum(diff(mz) < 0.1)
	if( z / x$point_count[k] < 0.05 ) {
		x$mz <- as.integer(round(x$mz))
		res <- .cdffix(x)
		return(res)
	}
	else {
		return(NULL)
	}
}

# max_assigned: number of m/z values assigned to the same integer (to detect
#              high mass accuracy)
.cdffix <- function(x, max_assigned=3L) {
	.Call("cdffix", x, max_assigned, PACKAGE="TargetSearch")
}

# function to open a CDF file and perform all checks:
#  1. m/z must be integer
#  2. intensity must be integer
#  3. data inconsistencies

.open.ncdf <- function(cdfFile, force=FALSE) {
	nc <- nc_open(cdfFile)
	ncData <- list()
	ncData$point_count <- ncvar_get(nc, "point_count")
	ncData$scanindex   <- ncvar_get(nc, "scan_index")
	ncData$intensity   <- ncvar_get(nc, "intensity_values")
	ncData$mz          <- ncvar_get(nc, "mass_values")
	ncData$rt          <- ncvar_get(nc, "scan_acquisition_time")
	nc_close(nc)

	if(any(ncData$scanindex < 0)) {
		message('Error:')
		message(sprintf('The NetCDF file %s seems to be corrupted.', cdfFile))
		stop('Unable to processs file. Aborting.')
	}

	ncData$intensity <- as.integer(ncData$intensity)

	ncData <- .check.mz.precision(ncData, force)
	if(is.null(ncData)) {
		stop(paste("Error processing file '", cdfFile, "'. It seems to contains",
			"high mass accuracy data.", sep=" "))
	}

	return(ncData)
}

# vim: set ts=4 sw=4:
