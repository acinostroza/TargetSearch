\name{TSExample}
\alias{TSExample}
\alias{RImatrix}
\alias{corRI}
\alias{metabProfile}
\alias{peakData}
\alias{refLibrary}
\alias{sampleDescription}
\docType{data}
\title{Example GC-MS data for TargetSearch Package}
\description{
    A \code{TargetSearch} example GC-MS data. This datasets contains TargetSearch object
    examples generated from a E.coli salt stress experiment (See package TargetSearchData).
}
\usage{
    data(TSExample)
}
\format{
    The data contains the following objects:
    \describe{
        \item{sampleDescription}{a \code{tsSample} object. The sample description.}
        \item{refLibrary}{a \code{tsLib} object. The reference library.}
        \item{rimLimits}{a \code{tsRim} object. The RI markers definition.}
        \item{RImatrix}{a matrix object. The retention time of the RI markers.}
        \item{corRI}{a matrix object. The sample RI.}
        \item{peakData}{a \code{tsMSdata} object. The intensities and RIs of all the
            masses that were searched for.}
        \item{metabProfile}{a \code{tsProfile} object. The metabolite profile.}
    }
}
\details{
    This dataset contain only the objects. The actual source files are provided
    by the package TargetSearchData.
}

\seealso{
    \code{\link[TargetSearch]{ImportLibrary}},
    \code{\link[TargetSearch]{ImportSamples}},
    \code{\link[TargetSearch]{ImportFameSettings}},
}
\examples{
require(TargetSearchData)

## The directory with the NetCDF GC-MS files
cdfpath <- file.path(find.package("TargetSearchData"), "gc-ms-data")
cdfpath
list.files(cdfpath)
samp.file <- file.path(cdfpath, "samples.txt")
rim.file  <- file.path(cdfpath, "rimLimits.txt")
lib.file  <- file.path(cdfpath, "library.txt")

# import files from package
sampleDescription <- ImportSamples(samp.file, CDFpath = cdfpath, RIpath = ".")
refLibrary        <- ImportLibrary(lib.file)
rimLimits         <- ImportFameSettings(rim.file, mass = 87)
# perform RI correction
RImatrix          <- RIcorrect(sampleDescription, rimLimits, massRange = c(85,320),
                   IntThreshold = 25, pp.method = "ppc", Window = 15)
# update median RI
refLibrary        <- medianRILib(sampleDescription, refLibrary)
# get the sample RI
corRI             <- sampleRI(sampleDescription, refLibrary, r_thres = 0.95)
# obtain the peak Intensities of all the masses in the library
peakData          <- peakFind(sampleDescription, refLibrary, corRI)
# make a profile of the metabolite data
metabProfile      <- Profile(sampleDescription, refLibrary, peakData, r_thres = 0.95)

# show the metabolite profile
profileInfo(metabProfile)
# show the matrix intensities
Intensity(metabProfile)
}
\keyword{datasets}
