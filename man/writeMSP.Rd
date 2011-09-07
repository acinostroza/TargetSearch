\name{writeMSP}
\alias{writeMSP}
\title{Save spectra in MSP format to be visualized in NIST}
\description{
  This function creates MSP format file from peak intensities that
  can be viewed with NIST.
}
\usage{
writeMSP(metlib, metprof, file, append = FALSE)
}
\arguments{
  \item{metlib}{ A \code{tsLib} object. A metabolite library. }
  \item{metprof}{ A \code{tsProfile} object. Usually the ouput of \code{\link{Profile}} or
	 \code{\link{ProfileCleanUp}} functions. }
  \item{file}{ A string naming the output file. }
  \item{append}{ Logical. If \code{TRUE} the results will be appended to \code{file}.
  Otherwise, it will overwrite the contents of \code{file}.}
}
\author{Alvaro Cuadros-Inostroza}
\seealso{ \code{\link{peakFind}}, \code{\link{Profile}}, \code{\link{ProfileCleanUp}},
\code{\linkS4class{tsLib}}, \code{\linkS4class{tsMSdata}}, \code{\linkS4class{tsProfile}} }
