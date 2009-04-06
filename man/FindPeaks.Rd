\name{FindPeaks}
\alias{FindPeaks}
\title{ Extract peaks from chromatogram files}
\description{
  This function extracts the maximum intensity of a list of masses in a given
  RI window.
}
\usage{
FindPeaks(my.files, refLib, columns = c("SPECTRUM", "RETENTION_TIME_INDEX"),
          showProgressBar = FALSE)
}
\arguments{
  \item{my.files}{ A character vector naming files to be searched. }
  \item{refLib}{ A numeric matrix with three columns. The second column contains
    the masses and the first and third column contains the RI limits. }
  \item{columns}{ A numeric vector with the positions of the columns \code{SPECTRUM} and
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

my.files <- RIfiles(sampleDescription)
# make a three column matrix: lower RI, mass, upper RI
refLib   <- refLib(refLibrary)
head(refLib)

# extract the peaks
peaks <- FindPeaks(my.files, refLib)
    
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\note{ This is an internal function not intended to be invoked directly. }
\seealso{ \code{\link{medianRILib}}, \code{\link{sampleRI}}, 
\code{\link{peakFind}},  \code{\linkS4class{tsMSdata}} }
