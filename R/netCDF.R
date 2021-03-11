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
#'
#' Convert from NetCDF format 3 into a custom TargetSearch NetCDF format 4. The
#' new NetCDF just contains the raw data in a matrix format in order to allow
#' easier and faster data manipulation.
#'
#' Starting from version 1.42.0, TargetSearch introduces a custom NetCDF file which
#' is used for faster and easier data manipulation. This means, ion traces within a
#' retention time can be quickly extracted, which if often required before
#' plotting. Formerly, this process required parsing the whole file before the data
#' could be extracted.
#'
#' The function only takes one file at the time, to convert many files at the
#' same time see the function [ncdf4_convert_from_path()] or the high level method
#' [ncdf4Convert()]. Alternativelly, you can call this function in a loop or
#' using the `lapply` family of functions.
#'
#' Keep in mind this function is intended for internal use (or advanced users); it is
#' exported for convenience. Using the method [ncdf4Convert()] is recommended.
#'
#' @section File structure:
#' The structure of the NetCDF format 4 is straightforward and the variables and
#' attributes are self-evident. The following variables are defined.
#'
#' * `retention_time` is a vector representing the retention time in seconds (double).
#' * `retention_index` is a vector representing the retention time indices (double).
#'   If missing, then the variable contains zeros. Its length is equal to the length
#'   of `retention_time`.
#' * `mass_range` is vector of length two containing the minimum and maximum m/z
#'   values (integer).
#' * `intensity` is matrix of intensity values (integer) where columns represent
#'   ion traces and rows are scans. The dimensions are length of "retention time"
#'   times the number of ions, ie, mass max - mass min + 1.
#'
#' In addition, the following attributes are defined. Note that only `creator` and
#' `version` are mandatory.
#'
#' * `creator` a string equal to "TargetSearch" (for indentification purposes).
#' * `version` file format version (string). Currently "1.1".
#' * `time_corrected` a flag (short integer) to indicate RI correction.
#' * `baseline_corrected` a flag (short integer) to indicate that the file was
#'   baseline corrected by TargetSearch.
#'
#' @param cdfFile The NetCDF file to be converted
#' @param outFile The new output file. If \code{NULL}, it replaces the \code{cdfFile}'s
#'   file extension (which should be \code{.cdf}) by \code{.nc4}. If the file
#'   extension is not \code{.cdf}, then \code{.nc4} is just appended.
#' @param force Logical. Set to \code{TRUE} to allow file overwrites, for example
#'   if the destination file still exists, in which case a warning is thrown. Default to \code{FALSE}.
#' @param baseline Logical. Whether or not baseline correct the input file.
#' @param \dots extra options passed to [baseline()].
#'
#' @note
#' Currently, it is not possible to reconstruct the original NetCDF file from the
#' converted file, especially if nominal mass or baseline correction was applied.
#' On the other hand, if the NetCDF files are exported from custom chromatogram
#' files (such as thermo raw files or LECO peg files), then the NetCDF 3 files
#' can be deleted safely as there is always a way to recover them.
#'
#' @return A string. The path to the converted file or invisible.
#' @seealso [ncdf4Convert()], [ncdf4_convert_from_path()], [baseline()]
#' @export
#' @md
#' @author Alvaro Cuadros-Inostroza
#' @examples
#' require(TargetSearchData)
#'
#' # get files from package TargetSearchData
#' cdfpath <- file.path(find.package("TargetSearchData"), "gc-ms-data")
#'
#' # choose any file
#' cdf <- file.path(cdfpath, '7235eg04.cdf')
#' nc4 <- '7235eg04.nc4' # save file in current path
#'
#' # run the function
#' ret <- ncdf4_convert(cdf, nc4)
#'
#' # the output should match the output file
#' stopifnot(ret == nc4)
#'
#' # Use mapply to convert many files at the same time.
#' cdf <- paste0('7235eg0', 6:8, '.cdf')
#' nc4 <- paste0('7235eg0', 6:8, '.nc4')
#' ret <- mapply(ncdf4_convert, file.path(cdfpath, cdf), nc4)
#' stopifnot(ret == nc4)
#'
`ncdf4_convert` <-
function(cdfFile, outFile=NULL, force=FALSE, baseline=FALSE, ...)
{
	# parameter assertions
	assert_that(is.string(cdfFile))
	assert_that(is.null(outFile) || is.string(outFile))
	assert_that(is.flag(force))
	assert_that(is.flag(baseline))

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
#' retention markers found by [RIcorrect()], or to force the markers time
#' if, for example, the RI markers were not co-injected with the biological
#' samples. It wraps around [rt2ri()]
#'
#' This function is similar to [fixRI()], with the difference that is acts
#' upon a single file, whereas [fixRI()] requires a [tsSample-class]
#' object.
#'
#' @md
#' @param cdfFile Path to the CDF file
#' @param observed The observed RI markers retention times'.
#' @param standard The RI of said markers.
#' @return Returns `invisible`
#' @note
#' This function is meant to be used internally. It is exposed for convenience.
#' @seealso [fixRI()], [RIcorrect()]
#' @author Alvaro Cuadros-Inostroza
#'
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

#' Convert CDF files to CDF4 from a path automatically
#'
#' Convert from NetCDF format 3 into a custom TargetSearch NetCDF format 4
#' automatically by scanning for CDF-3 files in given path and calling
#' the function [ncdf4_convert()] on them.
#'
#' This function simply wraps around [ncdf4_convert()]. It searches for
#' CDF-3 files and converts them to CDF-4. Note that the search is not
#' recursive.
#'
#' By default, the function scans the current path and saves them in
#' the same path if `out_path` is not specified.
#'
#' @param cdf_path the input path to scan for
#' @param out_path the output path in which the files will be saved
#' @param \dots extra options passed to [ncdf4_convert()], which in
#'    can be passed to [baseline()]
#'
#' @return a character vector of generated files or invisible.
#' @seealso [ncdf4_convert()], [baseline()]
#' @export
#' @author Alvaro Cuadros-Inostroza
#' @md
#' @examples
#' \dontrun{
#' # get files from package TargetSearchData
#' require(TargetSearchData)
#' cdfpath <- file.path(find.package("TargetSearchData"), "gc-ms-data")
#' ncdf4_convert_from_path(cdfpath, ".")
#' }
#'
`ncdf4_convert_from_path` <-
function(cdf_path=".", out_path=cdf_path, ...)
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

################################################################################

#' Extract data ranges from a NetCDF file format 4
#'
#' A flexible and convenient function to extract raw data from a NetCDF file format
#' 4 using time ranges and m/z ranges or values. This is a better (and faster)
#' alternative to the old [peakCDFextraction()] function, which reads the whole CDF
#' into memory, specially when only sections of the CDF are needed.
#'
#' The function takes a NetCDF format 4 generated by "TargetSearch" and extracts
#' raw intensity values from given m/z ion traces within a given time range. The
#' time range can be in seconds or arbitrary retention time index units. For
#' the latter case, the function expects a time corrected CDF file.
#'
#' If the given time range is out of range, a `NULL` value will be returned. In contrast,
#' if the m/z values are out of range, then zeros will be returned for out of range masses
#' (provided that the time range is not out of range). If `timeRange` is missing, then
#' the whole time range. Similarly, if `massRange` is missing, then all the masses are
#' extracted. If both are missing, the function behaves as [peakCDFextraction()], but
#' the output (a named list) uses slightly different names.
#'
#' The NetCDF must be have been previously converted to the custom "TargetSearch" format,
#' otherwise an error will be raised. See [ncdf4_convert()] for the convertion.
#'
#' @param cdfFile A path to a NetCDF file format 4.
#' @param massValues A numeric vector representing m/z values.
#' @param timeRange A numeric vector of length 2 representing the lower and upper
#'   time limits.
#' @param useRT Logical. If `TRUE`, the time range is in seconds, otherwise if `FALSE`
#'   (default), the time range is inretention time index units (`TRUE`).
#'
#' @return
#' A named list with the following components.
#'
#' * `Time` Numeric vector: the retention time in seconds
#' * `Index` Numeric vector: the retention time index (or zero if the file was
#'          was not retention time corrected
#' * `Intensity` Matrix: Rows are the retention times (or scans) and columns are masses.
#' * `massRange` Numeric vector of length 2: the mass range of the CDF
#'
#' @note
#' An error will be produced for invalid files. Also, this function is intended
#' to be used internally, but it is exposed for convinience, for example, to
#' create custom plots.
#'
#' @seealso [ncdf4_convert()], [peakCDFextraction()]
#' @export
#' @md
#' @author Alvaro Cuadros-Inostroza
#' @examples
#' \dontrun{
#' # set a NCDF-4 file
#' nc4file <- "/path/to/netcdf4.nc4"
#'
#' # extract all data (behaves like peakCDFextraction)
#' data <- ncdf4_data_extract(nc4file)
#'
#' # extract only certain m/z values
#' data <- ncdf4_data_extract(nc4file, massValues=c(116, 192))
#'
#' # to use mass ranges, use the colon (:) operator for example
#' data <- ncdf4_data_extract(nc4file, massValues=c(120:130, 200, 203:209))
#'
#' # restrict extraction to a retention index interval
#' data <- ncdf4_data_extract(nc4file, massValues=c(116, 192),
#'                            timeRange=c(200000, 220000))
#'
#' # same, but using retention time in seconds.
#' data <- ncdf4_data_extract(nc4file, massValues=c(116, 192),
#'                            timeRange=c(200, 300), useRT=TRUE)
#' }
`ncdf4_data_extract` <-
function(cdfFile, massValues, timeRange, useRT=FALSE)
{
    where <- function(x, k) which(x > k[1] & x < k[2])
    assert_that(is.string(cdfFile))
    nc_info <- .ncdf_info(cdfFile)

    if(nc_info$creator != 'TargetSearch')
        stop("Invalid CDF format")
    if(nc_info$time_corrected == 0 & useRT == FALSE)
        stop("Unable to extract RAW data in a non-time-corrected file")

    nc <- nc_open(cdfFile)
    on.exit(nc_close(nc))
    RT <- ncvar_get(nc, 'retention_time')
    RI <- ncvar_get(nc, 'retention_index')
    MR <- ncvar_get(nc, 'mass_range')

    if(missing(timeRange)) {
        index <- seq(RT)
    } else {
        assert_that(is.flag(useRT))
        assert_that(is.numeric(timeRange), length(timeRange) == 2)
        index <- if(useRT) where(RT, timeRange) else where(RI, timeRange)
    }

    if(length(index) == 0)
        return(NULL)

    if(missing(massValues)) {
        z <- ncvar_get(nc, 'intensity')[index,,drop=FALSE]
        massValues <- seq(MR[1], MR[2])
    } else {
        assert_that(is.numeric(massValues))
        massValues <- as.integer(massValues)

        z <- sapply(massValues, function(x) {
            if(x < MR[1] | x > MR[2])
                numeric(length(index))
            else
                ncvar_get(nc, 'intensity', start=c(index[1], x - MR[1] + 1), count=c(length(index), 1))
        })

        if(!is.matrix(z))
            z <- matrix(z, length(index),  length(massValues))
    }
    colnames(z) <- as.character(massValues)
    list(Time=RT[index], Index=RI[index], Intensity=z, massRange=MR)
}

#' plot peak from NetCDF-4 files
#'
#' @param obj list of files of tsSample object
#' @param massValues the m/z values (two-element vector)
#' @param timeRange the time range (two-element vector)
#' @param useRT use RT as time units (RI if FALSE)
#' @param showRT use RT as y-axis (RI if FALSE)
#' @param plot should a plot be drawn?
#' @param ... extra plotting parameters
#' @return a matrix with the data used for plotting
`ncdf4_plot_peak` <-
function(obj, massValues, timeRange, useRT=TRUE, showRT=useRT, plot=TRUE, ...)
{
    if(inherits(obj, 'tsSample')) {
        cdf <- CDFfiles(obj)
        names(cdf) <- sampleNames(obj)
    } else if(is.character(obj)) {
        cdf <- obj
    } else {
        stop('Object must inherit `tsSample` or be a character vector')
    }
    assert_that(is.flag(useRT))
    assert_that(is.flag(showRT))
    assert_that(is.flag(plot))

    optfun <- function(j, ...) {
        J <- function(x, j) x[ ((j - 1) %% length(x)) + 1]
        optnames <- c('scol', 'stype', 'slty', 'slwd', 'spch', 'sbg', 'scex')
        opts <- list(...)
        if(j == 0) {
            opts <- opts[ ! names(opts) %in% optnames ]
            return( opts[!duplicated(names(opts))] )
        }
        for(name in optnames) {
            if(!is.null(opts[[ name ]])) {
                opts[[ substring(name, 2)]] <- J(opts[[ name ]], j)
                opts[[name]] <- NULL
            }
        }
        opts[!duplicated(names(opts))]
    }

    dat <- lapply(cdf, ncdf4_data_extract, massValues, timeRange, useRT)
    if(is.null(names(dat)))
        names(dat) <- basename(cdf)

    null <- sapply(dat, is.null)
    if(all(null)) {
        warning("Nothing to plot")
        return(invisible())
    }

    if(!plot)
        return(invisible(dat))

    unit <- if(showRT) 'Time' else 'Index'
    xlim <- range(sapply(dat[!null], function(x) range(x[[ unit ]])))
    ylim <- range(sapply(dat[!null], function(x) range(x[[ 'Intensity' ]])))
    M <- length(massValues)

    opt <- optfun(0, x=1, type='n', ..., xlim=xlim, ylim=ylim, xlab=unit, ylab='Intensity')
    do.call('plot', opt)

    for(j in seq_along(dat)) {
        x <- dat[[j]]
        if(is.null(x))
            next
        opt <- optfun(j, x=x[[ unit ]], y=x[[ 'Intensity' ]], ...)
        do.call('matlines', opt)
    }

    return(invisible(dat))
}

# vim: set ts=4 sw=4:
