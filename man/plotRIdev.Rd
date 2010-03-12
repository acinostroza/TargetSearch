\name{plotRIdev}
\alias{plotRIdev}
\alias{plotAllRIdev}
\title{ Plot Retention Time Index Deviation }
\description{
\code{plotRIdev} plots the Retention Time Index Deviation of a given set of metabolites.
\code{plotAllRIdev} saves the plots of the RI deviations of all the metabolites in the library object
into a PDF file.
}
\usage{
plotRIdev(Lib, peaks, libId = 1)

plotAllRIdev(Lib, peaks, pdfFile, width = 8, height = 8, ...)
}
\arguments{
  \item{Lib}{ A \code{tsLib} object created by \code{\link{ImportLibrary}} function. }
  \item{peaks}{ A \code{tsMSdata} object. See \code{\link{peakFind}}.}
  \item{libId}{A numeric vector providing the indices of the metabolites to plot.}
  \item{pdfFile}{A file name where the plot will be saved. Only \code{plotAllRIdev}. }
  \item{width, height}{The width and height of the plots in inches. Only \code{plotAllRIdev}. }
  \item{\dots}{Further options passed to \code{\link{pdf}}. }
}
\examples{
require(TargetSearchData)
data(TargetSearchData)

# get RI file path
RI.path <- file.path(.find.package("TargetSearchData"), "gc-ms-data")
# update RI file path
RIpath(sampleDescription) <- RI.path

peakData <- peakFind(sampleDescription, refLibrary, corRI)

# Plot RI deviation of metabolite "Valine"
grep("Valine", libName(refLibrary)) # answer: 3
plotRIdev(refLibrary, peakData, libId = 3)

# Plot an RI deviation overview of the first nine metabolites
plotRIdev(refLibrary, peakData, libId = 1:9)

# Save all RI deviation into a pdf file
plotAllRIdev(refLibrary, peakData, pdfFile = "RIdeviations.pdf")

}

\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{ImportLibrary}}, \code{\linkS4class{tsLib}},
\code{\linkS4class{tsMSdata}}, \code{\link{pdf}} }
\keyword{ hplot }

