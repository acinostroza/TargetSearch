\name{plotSpectra}
\alias{plotSpectra}
\alias{plotAllSpectra}
\title{ Plot a Spectra Comparison}
\description{
\code{plotSpectra} plots a contrast between the reference spectra and the median spectra of a
given metabolite in the library. \code{plotAllRIdev} saves the plots of the median-reference
spectra comparisons of all the metabolites in the reference library into a PDF file.
}
\usage{
plotSpectra(Lib, peaks, libId = 1, type = "ht")

plotAllSpectra(Lib, peaks, type = "ht", pdfFile, width = 8, height = 8, ...)
}
\arguments{
  \item{Lib}{ A \code{tsLib} object created by \code{\link{ImportLibrary}} function. }
  \item{peaks}{ A \code{tsMSdata} object. See \code{\link{peakFind}}.}
  \item{libId}{ A numeric vector providing the indices of the metabolites
	to plot. }
  \item{type}{ The type of the plot. Options are \code{"ht"}, head-tail plot, \code{"ss"},
	side by side plot, and \code{"diff"}, spectrum difference plot. }
  \item{pdfFile}{ A file name where the plot will be saved. Only \code{plotAllRIdev}. }
  \item{width, height}{ The width and height of the plots in inches. Only \code{plotAllRIdev}. }
  \item{\dots}{ Further options passed to \code{pdf}. }
}
\details{
    The median spectra is obtained by computing the median intensity of every ion across the
    samples. The median and the reference spectra values are scaled to vary between 0 and 999
    in order to make them comparable.
}
\examples{
require(TargetSearchData)
data(TargetSearchData)

# get RI file path
RI.path <- file.path(find.package("TargetSearchData"), "gc-ms-data")
# update RI file path
RIpath(sampleDescription) <- RI.path

peakData <- peakFind(sampleDescription, refLibrary, corRI)

# Plot a comparison RI deviation of metabolite "Valine"
grep("Valine", libName(refLibrary)) # answer: 3
plotSpectra(refLibrary, peakData, libId = 3, type = "ht")

# Plot the spectra "side by side"
plotSpectra(refLibrary, peakData, libId = 3, type = "ss")

# Plot the spectra difference
plotSpectra(refLibrary, peakData, libId = 3, type = "diff")
}

\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\linkS4class{tsLib}}, \code{\linkS4class{tsMSdata}},\code{\link{pdf}} }
\keyword{ hplot }

