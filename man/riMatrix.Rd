\name{riMatrix}
\alias{riMatrix}
\title{Retention Time Index Matrix}
\description{
  A function to search for retention index RI markers.
}
\usage{
riMatrix(samples, rim)
}
\arguments{
  \item{samples}{ A \code{tsSample} object created by \code{\link{ImportSamples}} function. }
  \item{rim}{ A \code{tsRim} object. See \code{\link{ImportFameSettings}}.}
}
\details{
    This function works similar to \code{\link{RIcorrect}}, but searches
    for RI markers in RI files (not in CDF files). Can be used to retrieve
    the retention times of RI markers in already processed files.

    Note that it does not perform any RI adjusment. See \code{\link{fixRI}}.
}
\note{
    In case the RI files are in text format and their column names are not standard (for
    example, when the files were generated with another software), use the global option
    \code{'TS_RI_columns'} or transform the RI files to \code{\link{TargetSearch}} binary
    format. See the documentation in function \code{\link{text2bin}}.
}
\value{
    A retention time matrix of the found retention time markers. Every column represents a
    sample and rows RT markers.
}
\examples{
require(TargetSearchData)
# import refLibrary, rimLimits and sampleDescription.
data(TSExample)
# get the CDF files
cdfpath <- file.path(find.package("TargetSearchData"), "gc-ms-data")

# select a subset of samples
smp <- sampleDescription[1:4]

# update the CDF path
CDFpath(smp) <- cdfpath

# make a copy of the RI markers object
rim <- rimLimits

# run RIcorrect
RImat <- RIcorrect(smp, rim, massRange = c(85,320), writeCDF4path=FALSE,
           Window = 15, pp.method = "ppc", IntThreshold = 50)

# extract the retention times of the markers
RImat2 <- riMatrix(smp, rim)

# both matrices should be equal
stopifnot( all.equal(RImat, RImat2, tolerance=1e-8) )
}
\author{Alvaro Cuadros-Inostroza}
\seealso{ \code{\link{RIcorrect}}, \code{\link{FAMEoutliers}},\code{\link{ImportSamples}},
 \code{\link{ImportFameSettings}}, \code{\link{fixRI}}}
