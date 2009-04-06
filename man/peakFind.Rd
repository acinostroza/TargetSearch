\name{peakFind}
\alias{peakFind}
\title{ Intensities and RI matrices }
\description{
  This function returns a list of the intensities and RI matrices that were searched.
}
\usage{
peakFind(samples, Lib, cor_RI, columns = c("SPECTRUM", "RETENTION_TIME_INDEX"),
         showProgressBar = FALSE)
}
\arguments{
  \item{samples}{ A \code{tsSample} object created by \code{ImportSamples} function. }
  \item{Lib}{ A \code{tsLib} object created by \code{ImportLibrary} function with
	 corrected RI values.	See \code{medianRILib}. }
  \item{cor_RI}{ A matrix of correlating selective masses RI for every sample.
	 See \code{sampleRI}. }
  \item{columns}{ A numeric vector with the column positions of \code{SPECTRUM} and
	\code{RETENTION_TIME_INDEX} or a character vector with the header names of those
	columns. }
	\item{showProgressBar}{Logical. Should the progress bar be displayed?}
}
\value{
	A \code{tsMSdata} object.
}
\examples{
require(TargetSearchData)
data(TargetSearchData)

# get RI file path
RI.path <- file.path(.find.package("TargetSearchData"), "gc-ms-data")
# update RI file path
RIpath(sampleDescription) <- RI.path

peakData <- peakFind(sampleDescription, refLibrary, corRI)
# show peak Intensities. 
head(Intensity(peakData))

# How to get intensities for a particular metabolite
#
# make a library index using top masses
libId <- libId(refLibrary, sel = FALSE)
# get the peak intensities of Metabolite 1, for example, of every mass
int.1 <- Intensity(peakData)[libId == 1,]
# this assigns the mass values to the row names of int.1
rownames(int.1) <- topMass(refLibrary)[[1]]

}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{ImportSamples}}, \code{\link{ImportLibrary}}, \code{\link{medianRILib}},
	 \code{\link{sampleRI}}, \code{\linkS4class{tsMSdata}}, \code{\linkS4class{tsLib}},
	 \code{\linkS4class{tsSample}} }