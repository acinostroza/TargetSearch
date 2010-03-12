\name{plotFAME}
\alias{plotFAME}
\title{ Plot a standard marker }
\description{
	Plots a given standard marker.
}
\usage{
plotFAME(samples, RImatrix, whichFAME)
}
\arguments{
  \item{samples}{ A \code{tsSample} object created by \code{ImportSamples} function. }
  \item{RImatrix}{ A retention time matrix of the found retention time markers. }
  \item{whichFAME}{ The retention marker to plot. Must be a number between \code{1}
	and the number of markers.}
}
\examples{
require(TargetSearchData)
data(TargetSearchData)
# plot Retention index standards 1 to 3
plotFAME(sampleDescription, RImatrix, 1)
plotFAME(sampleDescription, RImatrix, 2)
plotFAME(sampleDescription, RImatrix, 3)
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{RIcorrect}}, \code{\link{FAMEoutliers}}, 
 \code{\linkS4class{tsSample}} }
\keyword{ hplot }

