\name{peakCDFextraction}
\alias{peakCDFextraction}
\title{ NetCDF to R }
\description{
    This function reads a netcdf chromatogram (format 3 or 4) file and returns
    a list containing raw data (retention time, intensity matrix) for easy
    manipulation.
}
\usage{
peakCDFextraction(cdfFile, massRange)
}
\arguments{
  \item{cdfFile}{A character string naming a netcdf file.}
  \item{massRange}{Deprecated. The mass range is extracted automatically. }
}
\details{
    The \code{cdfFile} is a path to a CDF-3 or a CDF-4 format file. The CDF-4
    format file is a custom file used on \code{\link{TargetSearch}}.

    The \code{massRange} parameter is deprecated but kept for compatibility.
    The actual m/z range is detected automatically and returned to the user.
}
\value{
A list with the following components:
  \item{Time}{The retention time vector.}
  \item{Index}{The retention time index vector or NULL if unavailable.}
  \item{Peaks}{The intensity matrix. Rows are the retention times and columns are masses. The first
    column is the lower mass value and the last one is the higher mass. }
  \item{massRange}{The mass range.}
  \item{baselineCorrected}{Logical. It is \code{TRUE} if the file was baseline corrected by
      the function \code{\link{baseline}}. See note below.}
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\note{
    This function does not look for peaks, just extracts all the raw intensity values
    of the chromatogram file. Use \code{\link{NetCDFPeakFinding}} instead.

    The function does not know whether a CDF file was baseline corrected (for
    example by the GC-MS software vendor), unless it was performed by the
    function \code{\link{baseline}}. It is perfectly possible to have a
    baseline corrected file and the element \code{baselineCorrected} be \code{FALSE}.
}
\seealso{\code{\link{ncdf4_data_extract}} for an updated and more flexible version
   of this function, \code{\link{NetCDFPeakFinding}}, \code{\link{baseline}}.}
\examples{
### take a NCDF3 file from the package TargetSearchData
require(TargetSearchData)
# pick a CDF file and extract its data
cdffile <- tsd_cdffiles()[4]
peakData <- peakCDFextraction(cdffile)

# it also works for NCDF4 files (we need to convert the file first)
nc4file <- ncdf4_convert(cdffile, tempfile(fileext='.nc4'))
peakData2 <- peakCDFextraction(nc4file)

# clean-up
unlink(nc4file)
}
% vim: set ts=4 sw=4 et:
