# auxiliary functions to NetCDFPeakFinding and peakCDFextraction

# function to open a CDF file (version 3)
#  1. m/z can be integer or double. Data will be converted to nominal mass.
#  2. intensity must be integer. Check if an overflow could occurr (R uses 32-bit integers)
#  3. data inconsistencies

.open.ncdf <- function(cdfFile) {
	nc <- nc_open(cdfFile)
	ncData <- list()
	ncData$point_count <- ncvar_get(nc, "point_count")
	ncData$scanindex   <- ncvar_get(nc, "scan_index")
	ncData$intensity   <- ncvar_get(nc, "intensity_values")
	ncData$mz          <- ncvar_get(nc, "mass_values")
	ncData$rt          <- ncvar_get(nc, "scan_acquisition_time")
	nc_close(nc)

	if(any(ncData$scanindex < 0) | any(diff(ncData$scanindex) < 0)) {
		message('Error:')
		message(sprintf('The NetCDF file %s seems to be corrupted.', cdfFile))
		stop('Unable to processs file. Aborting.')
	}

	# check intensity range
	if(max(ncData$intensity) > 2147483647)
		stop(paste("File", cdfFile, ": Intensity values cannot be coerced as",
			"integers as they will create an integer overflow. Aborting."))

	# sanity checks
	ncData <- .ncdf_sanity(ncData)

	if(!is.integer(ncData$intensity))
		ncData$intensity <- as.integer(ncData$intensity)

	# check if we need to transform to nominal mass
	if(!is.integer(ncData$mz)) {
		ncData$mz <- round(ncData$mz)
		ncData <- .Call(c_nominal, ncData)
	}
	return(ncData)
}

###############################################################################
# functions to detect netCDF format 3 or 4

#' Extract meta info from a netCDF file
#'
#' Returns a list with the NetCDF file format, and extracts the 'creator',
#' 'version' and 'time_corrected' attributes.
#'
#' @param cdf Path to the NetCDF file
#'
#' @return a list with components \code{format} (the netCDF format),
#' \code{creator} ('Unknown' if undefined), \core{version} (file version or empty
#' string if undefined) and \code{time_corrected} (1, 0, or NA_integer_ if undefined)
#'
`.get_ncdf_info` <- function(cdf)
{
	nc <- nc_open(cdf)
	format   <- nc$format
	creator <- ncatt_get(nc, 0, 'creator')
	timecor <- ncatt_get(nc, 0, 'time_corrected')
	version <- ncatt_get(nc, 0, 'version')
	nc_close(nc)

	creator <- if(creator$hasatt) creator$value else "Unknown"
	version <- if(version$hasatt) version$value else ""
	timecor <- if(timecor$hasatt) timecor$value else NA_integer_
	list(format=format, creator=creator, version=version, time_corrected=timecor)
}

#' Checks that a cdfFile was created by TS
`.is_ts_ncdf4` <- function(cdfFile)
{
	nfo <- .get_ncdf_info(cdfFile)
	nfo$creator == 'TargetSearch'
}

#' Convert a NetCDF file format 3 to format 4
#
#' Convert a NetCDF format 3 into a custom TargetSearch NetCDF format 4.
#' The new NetCDF just contains an intensity matrix (time x m/z) in order
#' to allow easier and faster data manipulation.
#'
#' @param cdfFile The NetCDF file to be converted
#' @param outFile The new output file. If \code{NULL}, it replaces the \code{cdfFile}'s
#'   file extension by \code{.nc4}. Valid extensions are \code{.cdf} or \code{.nc}. If
#'   the file doesn't have a is valid extension, then \code{.nc4} is just appended.
#' @param massRange The \code{m/z} range. It is actually ignored but kept for compatibility
#' @param force Logical. Set to \code{TRUE} to allow overwrites. Default to \code{FALSE}
#'
#' @note
#' The generated CDF file is non-standard and very likely cannot be
#' used outside targetSearch. For instance cannot be used in AMDIS.
#' It is not possible reconstruct the original NetCDF file. On the other hand,
#' if the NetCDF files are exported from the custom vendor files, then
#' the NetCDF 3 files can be deleted safely (as long you keep your original
#' files).
#'
#' @return A string. The path to the converted file or invisible.
`convert_to_ncdf4` <-
function(cdfFile, outFile=NULL, massRange=NULL, force=FALSE)
{
	# first check that the file is not TS
	if(.is_ts_ncdf4(cdfFile)) {
		warning('File "', cdfFile, '" has already been converted')
		return(invisible(outFile))
	}

	if(is.null(outFile)) {
		outFile <- sprintf("%s.nc4", sub("\\.cdf$", "", cdfFile, ignore.case=TRUE))
	}

	if(cdfFile == outFile)
		stop('Intput and output files are the same. aborting...')

	if(file.exists(outFile) & !force) {
		warning('File `', outFile, "' exists. Set `force' to overwrite")
		return(invisible(outFile))
    }
	peaks <- peakCDFextraction(cdfFile, massRange)
	save_cdf4(outFile, peaks)

	invisible(outFile)
}

`.save_cdf4_internal` <-
function(cdf, retTime, Peaks, massRange, retIndex=NULL, chunksizes=c(500, 5))
{
	n  <- ncol(Peaks)
	if(length(retTime) != nrow(Peaks))
		stop('length of retTime and nrow of Peaks must be equal')

	# define dimensions
	time_dim  <- ncdim_def('time', '', seq_along(retTime), create_dimvar=FALSE)
	mass_dim  <- ncdim_def('mass', '', seq_len(ncol(Peaks)), create_dimvar=FALSE)
	range_dim <- ncdim_def('range', '', 1:2, create_dimvar=FALSE)

	# define variables
	int_var <- ncvar_def('intensity', 'count', list(time_dim, mass_dim), prec='integer', compression=1, chunksizes=chunksizes)
	RT_var  <- ncvar_def('retention_time', 'second', time_dim, prec='double', compression=1)
	RI_var  <- ncvar_def('retention_index', 'unit',  time_dim, prec='double', compression=1)
	mr_var  <- ncvar_def('mass_range', 'mz', range_dim, prec='integer')

	ncnew <- nc_create(cdf, list(mr_var, int_var, RI_var, RT_var), force_v4=TRUE)
	on.exit( nc_close(ncnew) )

	ncvar_put( ncnew, mr_var, massRange)
	ncvar_put( ncnew, int_var, Peaks)
	ncvar_put( ncnew, RT_var, retTime)

	ncatt_put( ncnew, 0, 'creator', 'TargetSearch')
	ncatt_put( ncnew, 0, 'version', '1.0')

	if(is.null(retIndex)) {
		ncvar_put( ncnew, RI_var, numeric(length(retTime)))
			ncatt_put( ncnew, 0, 'time_corrected', 0, prec='short')
	} else {
		ncvar_put( ncnew, RI_var, retIndex)
		ncatt_put( ncnew, 0, 'time_corrected', 1, prec='short')
	}
}

# save peak structure to cdf4
`save_cdf4` <- function(cdf, peaks)
{
	if(!all(c('Time', 'Peaks', 'massRange') %in% names(peaks)))
		stop("Invalid list peaks. Missing names")
	.save_cdf4_internal(cdf, peaks$Time, peaks$Peaks, peaks$massRange, peaks$Index)
}

# Description
#  Update the retention index on the CDF files. The parameters observed
#  and standard are the same as used in the rt2ri function.

`update_retention_index_ncdf4` <-
function(cdfFile, observed, standard)
{
	# first check that the file is not TS
	if(!.is_ts_ncdf4(cdfFile)) {
		stop('File "', cdfFile, '" is not a recognized TS format')
	}

	nc     <- nc_open(cdfFile, write=TRUE)
	on.exit(nc_close(nc))
	rtTime <- ncvar_get(nc, 'retention_time')
	riTime <- rt2ri(rtTime, observed, standard)
	ncvar_put(nc, 'retention_index', riTime)
	ncatt_put(nc, 0, 'time_corrected', 1, prec='short')
	invisible()
}

#' scans for CDF files and converts them to CDF4
#'
#' @param cdf_path the input path to scan for
#' @param out_path the output path in which the files will be saved
`convert_cdf_from_path` <-
function(cdf_path, out_path=cdf_path)
{
	# scans CDF file
	in_files  <- dir(cdf_path, pattern="\\.cdf$", full.names=TRUE, ignore.case=TRUE)
	if(length(in_files) == 0)
		stop(sprintf("No CDF files detected in dir `%s`", cdf_path))
	out_files <- sub("\\.cdf$", ".nc4", basename(in_files), ignore.case=TRUE)
	out_files <- file.path(out_path, out_files)

	mapply(convert_to_ncdf4, in_files, out_files)
	invisible(out_files)
}

# CDF sanity check
.ncdf_sanity <- function(ncdf)
{
	si <- ncdf$scanindex
	pc <- ncdf$point_count
	ln <- length(ncdf$mz)
	ns <- length(si)

	if(ln != length(ncdf$intensity))
		stop("`mz` and `intensity` do no have equal lengths")

	if(any(ns != c(length(pc), length(ncdf$rt))))
		stop("invalid CDF file. data lengths are not equal")

	# check point count
	if(ln != sum(pc))
	{
		warning("Invalid `point_count` data. Trying `scanindex` instead...")
		pc <- diff(c(si,ln))
		if(ln != sum(pc))
			stop("Invalid CDF. `point_count` can be recovered")
		ncdf$point_count <- pc
	}
	else if(!all(pc == diff(c(si,ln)))) {
		warning("Invalid `scanindex` data. Trying `point_count` instead")
		ncdf$scanindex <- c(0, cumsum(p[-m]))
	}
	ncdf
}


# vim: set ts=4 sw=4:
