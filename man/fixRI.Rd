\name{fixRI}
\alias{fixRI}
\title{Fixing Retention Time Index Correction}
\description{
    This function can be used to correct the detected retention time index (RI)
    markers or to manually force their location to specific retention times if,
    for example, the RI markers were not co-injected with the biological samples.
}
\usage{
fixRI(samples, rimLimits, RImatrix=NULL, sampleNames=NULL)
}
\arguments{
  \item{samples}{ A \code{tsSample} object created by \code{\link{ImportSamples}} function. }
  \item{rimLimits}{ A \code{tsRim} object. See \code{\link{ImportFameSettings}}.}
  \item{RImatrix}{Optional. A retention time matrix of the found retention time markers that
  was obtained after running \code{\link{RIcorrect}}.}
  \item{sampleNames}{Optional. A character vector naming the samples that are to be RI corrected.}
}
\details{
    Sometimes the retention index limits (see \code{\link{ImportFameSettings}})
    are not set correctly and you will have to run the peak detection and RI
    correction function (\code{\link{RIcorrect}}) again, which may take a
    long time specially if there are many samples.

    Instead, a simple approach is to fix the RI limits and use this function
    to correct the generated RI files. Since these files are much smaller
    than CDF files (chromatograms), this runs much faster.

    Other possibility is that the time positions of one
    or more RI markers were wrongly detected because there was just simply no peak or
    the RI markers where not co-injected in some samples. You could manually
    force the locations of the RI markers. This case is discussed in the
    ``RICorrection'' vignette.

    The behavior of this function depends on whether the parameter \code{RImatrix}
    is \code{NULL} or not. If \code{NULL}, the RI markers will be searched again
    (using the settings of \code{rimLimits} parameters, which you should have already
    fixed) and the resulting values will be used to correct the RI files.
    If it is a numeric matrix, then these values will be used to correct the
    RI files. Note that this matrix dimmensions are exactly \emph{m} samples (rows)
    times \emph{n} (columns) RI markers.

    \code{sampleNames} controls which samples will be corrected. If \code{NULL}
    then all samples will be corrected. It could be character vector (sample
    names) or a numeric vector representing the sample indexes.
}
\value{
    Invisible \code{NULL}. It prints the corrected samples.
}
\examples{
require(TargetSearchData)
# import refLibrary, rimLimits and sampleDescription.
data(TargetSearchData)
# get the CDF files
cdfpath <- file.path(find.package("TargetSearchData"), "gc-ms-data")

# select a subset of samples
smp <- sampleDescription[1:4]

# update the CDF path
CDFpath(smp) <- cdfpath

# make a copy of the RI markers object
rim <- rimLimits

# mess up the limits of marker 3 (real value is 369 seconds app.)
rimLimits(rim)[3,] <- c(375, 400)

# run RIcorrect
RImat <- RIcorrect(smp, rim, massRange = c(85,320),
           Window = 15, pp.method = "ppc", IntThreshold = 50)

# fix the limits of marker 3
rimLimits(rim)[3,] <- c(360, 400)

# you could run again RIcorrect, but this is faster
fixRI(smp, rim)

# get the RI matrix
RImat <- riMatrix(smp, rim)

# compare the values with the real ones (previously stored in RImatrix)
stopifnot( all.equal(RImat, RImatrix[,1:4], tolerance=1e-8) )

# manual adjusment or RI markers for sample 3.
# Warning: this is just an example to illustrate how to use this function.
#          don't do this unless you know what you're doing.
RImat[,3] <- c(252, 311, 369)
fixRI(smp, rim, RImat, 3)
}
\author{Alvaro Cuadros-Inostroza}
\seealso{ \code{\link{RIcorrect}}, \code{\link{FAMEoutliers}},\code{\link{ImportSamples}},
 \code{\link{ImportFameSettings}}}
