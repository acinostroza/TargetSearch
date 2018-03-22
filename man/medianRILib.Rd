\name{medianRILib}
\alias{medianRILib}
\title{ Median RI library correction }
\description{
Return a \code{tsLib} object with the median RI of the selective masses across samples.
}
\usage{
medianRILib(samples, Lib, makeReport = FALSE, pdfFile = "medianLibRep.pdf",
	columns = c("SPECTRUM", "RETENTION_TIME_INDEX", "RETENTION_TIME"),
    showProgressBar = FALSE)
}

\arguments{
  \item{samples}{ A \code{tsSample} object created by \code{\link{ImportSamples}} function. }
  \item{Lib}{ A \code{tsLib} object created by \code{\link{ImportLibrary}} function. }
  \item{makeReport}{ Logical. If \code{TRUE} will report the RI deviations for every metabolite
	in the library. }
  \item{pdfFile}{ The file name where the report will be saved. }
  \item{columns}{ A numeric vector with the positions of the columns \code{SPECTRUM},
    \code{RETENTION_TIME_INDEX}, and \code{RETENTION_TIME} or a character vector
    with the header names of those columns.}
	\item{showProgressBar}{Logical. Should the progress bar be displayed?}
}
\value{
	A \code{tsLib} object. It will update the slot \code{med_RI} which contains
 the median RI of every	searched metabolite.
}
\examples{
require(TargetSearchData)
data(TSExample)

# get RI file path
RI.path <- file.path(find.package("TargetSearchData"), "gc-ms-data")
# update RI file path
RIpath(sampleDescription) <- RI.path
# Import Library
refLibrary        <- ImportLibrary(file.path(RI.path,'library.txt'))
# update median RI
refLibrary        <- medianRILib(sampleDescription, refLibrary)

# perhaps you need to adjust the library RI of one metabolite and the allowed time
# deviation (first time deviation window)
libRI(refLibrary)[5] <- 306500
RIdev(refLibrary)[5,1] <- 2000

refLibrary        <- medianRILib(sampleDescription, refLibrary)

}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{ImportSamples}}, \code{\link{ImportLibrary}}, \code{\link{tsLib-class}} }
