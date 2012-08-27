\name{TargetSearch-deprecated}
\alias{TargetSearch-deprecated}
\alias{fixRIcorrection}
\alias{fixRIcorrection-deprecated}
\title{Deprecated Functions in Package \pkg{TargetSearch}}
\description{
  These functions are provided for compatibility with older versions of
  \pkg{TargetSearch}, and be defunct as soon as the next bioconductor
  release.
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
  \code{\link{fixRIcorrection}} has been replaced by \code{\link{fixRI}} which
  provides the same functionality and more.
}
\seealso{
  \code{\link{Deprecated}}
}
\keyword{misc}

