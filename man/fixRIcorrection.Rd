\name{fixRIcorrection}
\alias{fixRIcorrection}
\title{Manual Retention Time Index Correction}
\description{
    This function can be used to manually correct the detected retention time index (RI)
    or to force their location to specific time.
}
\usage{
fixRIcorrection(samples, rimLimits, RImatrix, sampleNames)
}
\arguments{
  \item{samples}{ A \code{tsSample} object created by \code{\link{ImportSamples}} function. }
  \item{rimLimits}{ A \code{tsRim} object. See \code{\link{ImportFameSettings}}.}
  \item{RImatrix}{A retention time matrix of the found retention time markers that
  was obtained after running \code{\link{RIcorrect}}.}
  \item{sampleNames}{A character vector naming the samples that are to be RI corrected.}
}                                                                                    
\details{
    This function should not be needed, unless the time positions of one
    or more RI markers were wrongly detected, which almost never occurs.
    If that were to happen, check and adjust the detection time limits and run
    again \code{\link{RIcorrect}} before using this function.

    The objects \code{samples}, \code{rimLimits} and \code{RImatrix} must be the
    same as those that were used in a previous call of \code{\link{RIcorrect}}.
}
\author{Alvaro Cuadros-Inostroza}
\seealso{ \code{\link{RIcorrect}}, \code{\link{FAMEoutliers}},\code{\link{ImportSamples}},
 \code{\link{ImportFameSettings}}}
