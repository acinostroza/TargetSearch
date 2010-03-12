\name{writeMSP}
\alias{writeMSP}
\title{Save spectra in MSP format to be visualized in NIST}
\description{
  This function creates MSP format file from peak intensities that
  can be viewed with NIST.
}
\usage{
writeMSP(lib, peaks, file, append = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lib}{ A \code{tsLib} object. }
  \item{peaks}{ A \code{tsMSdata} object. }
  \item{file}{ A string naming the output file. }
  \item{append}{ Logical. If \code{TRUE} the results will be appended to \code{file}.
  Otherwise, it will overwrite the contents of \code{file}.}
}
\author{Alvaro Cuadros-Inostroza}
\seealso{ \code{\link{peakFind}}, \code{\link{Profile}}, \code{\link{ProfileCleanUp}},
\code{\linkS4class{tsLib}}, \code{\linkS4class{tsMSdata}}, \code{\linkS4class{tsProfile}} }
