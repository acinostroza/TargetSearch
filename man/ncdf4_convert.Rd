% Please see file R/netCDF.R
\name{ncdf4_convert}
\alias{ncdf4_convert}
\title{Convert from a NetCDF file format 3 to format 4}
\usage{
ncdf4_convert(cdfFile, outFile = NULL, force = FALSE,
  baseline = FALSE, ...)
}
\arguments{
\item{cdfFile}{The NetCDF file to be converted}

\item{outFile}{The new output file. If \code{NULL}, it replaces the \code{cdfFile}'s
file extension (which should be \code{.cdf}) by \code{.nc4}. If the file
extension is not \code{.cdf}, then \code{.nc4} is just appended. If the path to
the file does not exist, it will be created automatically.}

\item{force}{Logical. Set to \code{TRUE} to allow file overwrites, for example
if the destination file still exists, in which case a warning is thrown. Default to \code{FALSE}.}

\item{baseline}{Logical. Whether or not baseline correct the input file.}

\item{\dots}{extra options passed to \code{\link[=baseline]{baseline()}}.}
}
\value{
A string. The path to the converted file or invisible.
}
\description{
Convert from NetCDF format 3 into a custom TargetSearch NetCDF format 4. The
new NetCDF just contains the raw data in a matrix format in order to allow
easier and faster data manipulation.
}
\details{
Starting from version 1.42.0, TargetSearch introduces a custom NetCDF file which
is used for faster and easier data manipulation. This means, ion traces within a
retention time can be quickly extracted, which if often required before
plotting. Formerly, this process required parsing the whole file before the data
could be extracted.

Note that this function only takes one file at the time. To convert many files at the
same time, see the function \code{\link[=ncdf4_convert_from_path]{ncdf4_convert_from_path()}} or the high level method
\code{\link[=ncdf4Convert]{ncdf4Convert()}}. Alternatively, you can call this function in a loop or
using the \code{\link[base:lapply]{lapply}} family of functions.

Keep in mind this function is intended for internal use (or advanced users); it is
exported for convenience. Using the method \code{\link[=ncdf4Convert]{ncdf4Convert()}} is recommended.
}
\note{
Currently, it is not possible to reconstruct the original NetCDF file from the
converted file, especially if nominal mass or baseline correction was applied.
On the other hand, if the NetCDF files are exported from custom chromatogram
files (such as thermo raw files or LECO peg files), then the NetCDF 3 files
can be deleted safely as there is always a way to recover them.
}
\section{File structure}{

The structure of the NetCDF format 4 is straightforward and the variables and
attributes are self-evident. The following variables are defined.
\itemize{
\item \code{retention_time} is a vector representing the retention time in seconds (double).
\item \code{retention_index} is a vector representing the retention time indices (double).
If missing, then the variable contains zeros. Its length is equal to the length
of \code{retention_time}.
\item \code{mass_range} is vector of length two containing the minimum and maximum m/z
values (integer).
\item \code{intensity} is matrix of intensity values (integer) where columns represent
ion traces and rows are scans. The dimensions are length of "retention time"
times the number of ions, i.e., mass max - mass min + 1.
}

In addition, the following attributes are defined. Note that only \code{creator} and
\code{version} are mandatory.
\itemize{
\item \code{creator} a string equal to "TargetSearch" (for identification purposes).
\item \code{version} file format version (string). Currently "1.1".
\item \code{time_corrected} (optional) a flag (short integer) to indicate RI correction.
      If missing it defaults to \code{false}.
\item \code{baseline_corrected} (optional) a flag (short integer) to indicate that the file was
      baseline corrected by TargetSearch. If missing it defaults to \code{false}.
}
}

\examples{
require(TargetSearchData)

# get files from package TargetSearchData
cdfpath <- tsd_data_path()

# choose any file
cdf <- file.path(cdfpath, '7235eg04.cdf')
nc4 <- '7235eg04.nc4' # save file in current path

# run the function
ret <- ncdf4_convert(cdf, nc4)

# the output should match the output file
stopifnot(ret == nc4)

# Use mapply to convert many files at the same time.
cdf <- paste0('7235eg0', 6:8, '.cdf')
nc4 <- paste0('7235eg0', 6:8, '.nc4')
ret <- mapply(ncdf4_convert, file.path(cdfpath, cdf), nc4)
stopifnot(ret == nc4)

}
\seealso{
\code{\link[=ncdf4Convert]{ncdf4Convert()}}, \code{\link[=ncdf4_convert_from_path]{ncdf4_convert_from_path()}}, \code{\link[=baseline]{baseline()}}
}
\author{
Alvaro Cuadros-Inostroza
}
