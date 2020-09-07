\name{peakFind}
\alias{peakFind}
\title{ Intensities and RI matrices }
\description{
  This function returns a list of the intensities and RI matrices that were searched.
}
\usage{
peakFind(samples, Lib, cor_RI, columns = NULL, showProgressBar = FALSE)
}
\arguments{
  \item{samples}{ A \code{tsSample} object created by \code{ImportSamples} function. }
  \item{Lib}{ A \code{tsLib} object created by \code{ImportLibrary} function with
	 corrected RI values. See \code{medianRILib}. }
  \item{cor_RI}{ A matrix of correlating selective masses RI for every sample.
	 See \code{sampleRI}. }
  \item{columns}{Either \code{NULL}, a \code{character} vector, or an \code{integer} vector. In
    most cases, leave it as \code{NULL}. This parameter is used to configure the column names or
    positions of RI text files. See the documentation on the \code{\link{text2bin}} function for
    further details.}
	\item{showProgressBar}{Logical. Should the progress bar be displayed?}
}
\value{
	A \code{tsMSdata} object.
}
\examples{
require(TargetSearchData)
data(TSExample)

# get RI file path
RI.path <- file.path(find.package("TargetSearchData"), "gc-ms-data")
refLibrary <- ImportLibrary(file.path(RI.path,"library.txt"))

# update RI file path
RIpath(sampleDescription) <- RI.path

peakData <- peakFind(sampleDescription, refLibrary, corRI)
# show peak Intensities. 
head(Intensity(peakData), 2)

# How to get intensities for a particular metabolite
# just select the identifier. Here extract the intensities
# for the first metabolite in the library
IntMatrix <- Intensity(peakData)[[1]]

}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{ImportSamples}}, \code{\link{ImportLibrary}}, \code{\link{medianRILib}},
	 \code{\link{sampleRI}}, \code{\linkS4class{tsMSdata}}, \code{\linkS4class{tsLib}},
	 \code{\linkS4class{tsSample}} }
