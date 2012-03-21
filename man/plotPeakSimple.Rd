\name{plotPeakSimple}
\alias{plotPeakSimple}
\title{ Plot peaks - simple interface }
\description{
	Plot selected ions in a given time range.
}
\usage{
plotPeakSimple(rawpeaks, time.range, masses, cdfFile = NULL, useRI = FALSE,
	 rimTime = NULL, standard = NULL, massRange = c(85, 500), ...)
}
\arguments{
  \item{rawpeaks}{ A two component list containing the retention time and the intensity matrices.
	See \code{\link{peakCDFextraction}}. }
  \item{time.range}{ The time range to plot in retention time or retention time index units to plot. }
  \item{masses}{ A vector containing the ions or masses to plot. }
  \item{cdfFile}{ The name of a CDF file. If a file name is specified, the ions will be extracted
	from there instead of using \code{rawpeaks}. }
  \item{useRI}{Logical. Whether to use Retention Time Indices or not.}
  \item{rimTime}{A retention time matrix of the found retention time markers. It is only used
	when \code{useRI} is \code{TRUE}.}
  \item{standard}{A numeric vector with RI values of retention time markers. It is only used
	when \code{useRI} is \code{TRUE}.}
  \item{massRange}{ A two component numeric vector with the scan mass range to extract. }
  \item{\dots}{ Further options passed to \code{\link{matplot}}. }
}
\note{This function used to named 'plotPeak'. This function was completely rewritten
so we kept the old version and renamed it 'plotPeakSimple'.}
\examples{
require(TargetSearchData)
data(TargetSearchData)

# update CDF path
CDFpath(sampleDescription) <- file.path(.find.package("TargetSearchData"), "gc-ms-data")

# Plot the peak "Valine" for sample number 1
grep("Valine", libName(refLibrary)) # answer: 3
# select the first file
cdfFile  <- CDFfiles(sampleDescription)[1]

# select "Valine" top masses
top.masses <- topMass(refLibrary)[[3]]

# plot peak from the cdf file
plotPeakSimple(cdfFile = cdfFile, time.range = libRI(refLibrary)[3] + c(-2000,2000),
    masses = top.masses, useRI = TRUE, rimTime = RImatrix[,1],
    standard = rimStandard(rimLimits), massRange = c(85, 500))

# the same, but extracting the peaks into a list first. This may be better if
# you intend to loop through several peaks.
rawpeaks <- peakCDFextraction(cdfFile, massRange = c(85,500))
plotPeakSimple(rawpeaks, time.range = libRI(refLibrary)[3] + c(-2000,2000),
    masses = top.masses, useRI = TRUE, rimTime = RImatrix[,1],
    standard = rimStandard(rimLimits), massRange = c(85, 500))

}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{\code{\link{plotPeak}}, \code{\link{RIcorrect}}, \code{\linkS4class{tsMSdata}}, \code{\linkS4class{tsRim}},
\code{\link{peakCDFextraction}}, \code{\link{matplot}} }
\keyword{ hplot }

