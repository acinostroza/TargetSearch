\docType{methods}
\name{ncdf4Convert,tsSample-method}
\alias{ncdf4Convert,tsSample-method}
\alias{ncdf4Convert}
\title{ Method for converting CDF-3 files to CDF-4  }
\description{
    \code{ncdf4Convert} is a high level method used to convert from
    CDF-3 to \sQuote{TargetSearch} new CDF-4 custom file format.
}
\usage{
    ncdf4Convert(obj, path, ...)
}
\arguments{
    \item{obj}{A \code{\linkS4class{tsSample}} object}
    \item{path}{A character vector representing file path in which the newly
        created CDF-4 files will be stored. If missing, the files will be
        save in the same directory as the source CDF-3 files}
    \item{\dots}{Extra arguments passed to \code{\link{ncdf4_convert}}, such as
        baseline correction options.}
}
\details{
    This is a high level interface to \code{\link{ncdf4_convert}} which
    uses a \code{\linkS4class{tsSample}} as input. The advantage of using
    this function is that it updates the paths to the CDF-4 files for
    downstream analysis.

    The parameter \code{path} can be used to change the original paths,
    which may be necessary if those files are on a read-only file-system.
    If this parameter is missing, the same paths are used. Note that
    the path are re-cycled to match the length of the samples.

    The \dots can be used to pass extra arguments to \code{\link{ncdf4_convert}}
    and to \code{\link{baseline}}. Refer to those man pages for details.
}
\value{
    A new \code{\linkS4class{tsSample}} object with updated CDF-4 file paths.
    The files are converted in the background.
}
\examples{
require(TargetSearchData)
data(TSExample)
# get the CDF files
cdfpath <- file.path(find.package("TargetSearchData"), "gc-ms-data")

# update the CDF path
CDFpath(sampleDescription) <- cdfpath

# tranform the CDF (the files are copied in the current directoru)
newSamples <- ncdf4Convert(sampleDescription, path=".")
}
\keyword{methods}
\keyword{tsSample}
\keyword{netcdf}
\author{Alvaro Cuadros-Inostroza}
% vim: set ts=4 sw=4 et:
