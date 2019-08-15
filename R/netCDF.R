# auxiliary functions to NetCDFPeakFinding and peakCDFextraction

#' Function to open a CDF3 file
#'
#' Opens and NetCDF format 3, normally exported from a software instrument,
#' and returns a list with the CDF data.
#'
#' The mass values can be integers or doubles. Because TargetSearch does not
#' use exact mass, they will converted to nominal mass by summing up
#' the intensity values closest to the nearest integer.
#'
#' The intensity values will be coerced to integers. If it is not possible,
#' the function will throw an error.
#'
#' Some sanity checks and performed and in some cases the data can be fixed,
#' otherwise an error will be thrown
#'
#' @param cdfFile A path to a CDF format 3 file
#' @return A list representing a CDF structure (which is quite standard)
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
#' 'version', 'time_corrected' and 'baseline_corrected' attributes.
#'
#' @param cdf Path to the NetCDF file
#'
#' @return A list with components \code{format} (the netCDF format),
#' \code{creator} ('Unknown' if undefined), \core{version} (file version or empty
#' string if undefined), \code{time_corrected} (1, 0, or NA_integer_ if undefined)
#' and \code{baseline_corrected} (1 if yes, 0 if not or unknown).
`.ncdf_info` <- function(cdf)
{
	nc <- nc_open(cdf)
	format   <- nc$format
	creator <- ncatt_get(nc, 0, 'creator')
	timecor <- ncatt_get(nc, 0, 'time_corrected')
	version <- ncatt_get(nc, 0, 'version')
	basecor <- ncatt_get(nc, 0, 'baseline_corrected')
	nc_close(nc)

	creator <- if(creator$hasatt) creator$value else "Unknown"
	version <- if(version$hasatt) version$value else ""
	timecor <- if(timecor$hasatt) timecor$value else NA_integer_
	basecor <- if(basecor$hasatt) basecor$value else 0
	list(format=format, creator=creator, version=version,
		 time_corrected=timecor, baseline_corrected=basecor)
}

#' Checks that a cdfFile was created by TS
#'
#' @param cdfFile Path to a CDF file
#' @return \code{TRUE} if is a TargetSearch file, \code{FALSE} otherwise.
`.is_ts_ncdf4` <- function(cdfFile)
{
	nfo <- .ncdf_info(cdfFile)
	nfo$creator == 'TargetSearch'
}

#' Convert from a NetCDF file format 3 to format 4
#
#' Convert from NetCDF format 3 into a custom TargetSearch NetCDF format 4.
#' The new NetCDF just contains an intensity matrix (time x m/z) in order
#' to allow easier and faster data manipulation.
#'
#' @param cdfFile The NetCDF file to be converted
#' @param outFile The new output file. If \code{NULL}, it replaces the \code{cdfFile}'s
#'   file extension (which should be \code{.cdf}) by \code{.nc4}. If the file
#'   extension is not \code{.cdf}, then \code{.nc4} is just appended.
#' @param force Logical. Set to \code{TRUE} to allow overwrites. Default to \code{FALSE}.
#' @param baseline baseline correct the CDF file
#' @param \dots extra options passed to baseline
#'
#' @note
#' This function is intended for internal use (or advanced users); it is exported
#' for convenience.
#'
#' The generated CDF file is non-standard and very likely cannot be
#' used outside targetSearch. For instance cannot be used by AMDIS.
#' Moreover, the mass values are converted to nominal mass (if they are not),
#' meaning there is data loss.
#'
#' Currently, it is not possible to reconstruct the original NetCDF file
#' from the converted file. On the other hand, if the NetCDF files are exported
#' from the custom vendor files, then the NetCDF 3 files can be deleted safely
#' (as long you keep your original chromatograms).
#'
#' @return A string. The path to the converted file or invisible.
#' @example
#' require(TargetSearchData)
#' cdfpath <- file.path(find.package("TargetSearchData"), "gc-ms-data")
#' cdf <- file.path(cdfpath, '7235eg04.cdf')
#' nc4 <- '7235eg04.nc4' # save file in current path
#' ret <- ncdf4_convert(cdf, nc4)
#' stopifnot(ret == nc4)
`ncdf4_convert` <-
function(cdfFile, outFile=NULL, force=FALSE, baseline=FALSE, ...)
{
	# extract information of CDF file
	nfo <- .ncdf_info(cdfFile)
	peaks <- NULL

	if(nfo$creator == 'TargetSearch') {
		if(baseline) {
			if(nfo$baseline_corrected == 0) {
				peaks <- peakCDFextraction(cdfFile)
				peaks <- baseline(peaks, ...)
			}
			else {
				warning('File "', cdfFile, '" has already been baseline corrected converted')
				return(invisible(outFile))
			}
		}
		else {
			warning('File "', cdfFile, '" has already been converted')
			return(invisible(outFile))
		}
	}

	if(is.null(outFile)) {
		outFile <- sprintf("%s.nc4", .trim_file_ext(cdfFile, 'cdf'))
	}

	if(cdfFile == outFile & !force)
		stop('Input and output files are the same. Set `force` to overwrite.')

	if(file.exists(outFile) & !force) {
		warning('File `', outFile, "' exists. Set `force' to overwrite")
		return(invisible(outFile))
    }
	if(is.null(peaks)) {
		peaks <- peakCDFextraction(cdfFile) # CDF3
		if(baseline)
			peaks <- baseline(peaks, ...)
	}

	ncdf4_write(outFile, peaks)
	invisible(outFile)
}

#' Internal function for saving data into a cdf4
#'
#' @param cdf Path to the CDF file.
#' @param retTime vector of retention times.
#' @param Peaks matrix of intensities.
#' @param massRange the mass range. A vector of length 2 (mz min, mz max).
#' @param retIndex the retention indices or NULL if not available.
#' @param baseline whether the intensity data was baseline corrected by TS.
#' @param chunksizes vector of chunk sizes for compression (see ncdf4)
`.ncdf4_write` <-
function(cdf, retTime, Peaks, massRange, retIndex=NULL, baseline=FALSE, chunksizes=c(500, 5))
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
	ncatt_put( ncnew, 0, 'version', '1.1')
	ncatt_put( ncnew, 0, 'baseline_corrected', as.integer(baseline[1]), prec='short')

	if(is.null(retIndex)) {
		ncvar_put( ncnew, RI_var, numeric(length(retTime)))
		ncatt_put( ncnew, 0, 'time_corrected', 0, prec='short')
	} else {
		ncvar_put( ncnew, RI_var, retIndex)
		ncatt_put( ncnew, 0, 'time_corrected', 1, prec='short')
	}
}

#' Save peak structure to cdf4
#'
#' Save a NCDF list to a CDF4 for. The structure is usually generated
#' by the function [peakCDFextraction()].
#' @param cdf A path to a CDF file to be written
#' @param peaks A list representing a NCDF4 structure. The list is generated
#'        by [peakCDFextraction()].
#' @note
#' This function is meant to be used internally. It is exposed for convenience.
`ncdf4_write` <- function(cdf, peaks)
{
	if(!all(c('Time', 'Peaks', 'massRange') %in% names(peaks)))
		stop("Invalid list peaks. Missing names")
	if(is.null(peaks$baselineCorrected))
		peaks$baselineCorrected <- FALSE
	.ncdf4_write(cdf, peaks$Time, peaks$Peaks, peaks$massRange, peaks$Index,
		peaks$baselineCorrected)
}

#' Update retention time index on a NCDF4 file
#'
#' Performs retention time index (RI) correction on a CDF file, using the
#' using the retention markers found by [RIcorrect()]. It wraps around
#' the function rt2ri()].
#'
#' @param cdfFile Path to the CDF file
#' @param observed The observed RI markers retention times'.
#' @param standard The RI of said markers.
#' @note
#' This function is meant to be used internally. It is exposed for convenience.
`ncdf4_update_ri` <-
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
#' @param ... extra options passed to ncdf4_convert
`ncdf4_convert_from_path` <-
function(cdf_path, out_path=cdf_path, ...)
{
	# scans CDF file
	in_files  <- dir(cdf_path, pattern="\\.cdf$", full.names=TRUE, ignore.case=TRUE)
	if(length(in_files) == 0)
		stop(sprintf("No CDF files detected in dir `%s`", cdf_path))
	out_files <- sub("\\.cdf$", ".nc4", basename(in_files), ignore.case=TRUE)
	out_files <- file.path(out_path, out_files)

	mapply(ncdf4_convert, in_files, out_files, MoreArgs=list(...))
	invisible(out_files)
}

# CDF sanity check
# checks that the point count match the respective scan indices. If not
# then uses one of those to compute the other. The checks are not exhaustive.
# A C-based check will be better.
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
		ncdf$scanindex <- c(0, cumsum(pc[-ns]))
	}
	ncdf
}

# vim: set ts=4 sw=4:
