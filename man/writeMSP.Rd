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
  \item{metprof}{ A \code{tsProfile} object. Usually the output of \code{\link{Profile}} or
	 \code{\link{ProfileCleanUp}} functions. }
  \item{file}{ A string naming the output file. }
  \item{append}{ Logical. If \code{TRUE} the results will be appended to \code{file}.
  Otherwise, it will overwrite the contents of \code{file}.}
}
\details{
  The generated file can be used to search in other library spectra using the
  software \emph{NIST Search Software}. This software is available here (only
  MS Windows is supported)
  \url{https://chemdata.nist.gov/dokuwiki/doku.php?id=chemdata:nistlibs}

  Note: the output file usually has extension \code{.msp}, which conflicts with
  MS Windows Software Patch files. Care is needed to not accidentally double
  click the file.
}
\author{Alvaro Cuadros-Inostroza}
\seealso{ \code{\link{peakFind}}, \code{\link{Profile}}, \code{\link{ProfileCleanUp}},
\code{\linkS4class{tsLib}}, \code{\linkS4class{tsMSdata}}, \code{\linkS4class{tsProfile}} }
\examples{
# here we use precomputed objects, in particular the reference library
#(refLibrary) and the metabolite profile (metabProfile).
data(TSExample)

# to call the function, simply use the aforementioned objects as inputs
# the file 'output_file.msp' can be open in NIST Search Software
writeMSP(refLibrary, metabProfile, file="output_file.msp")
}
