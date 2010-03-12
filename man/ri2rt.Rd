\name{ri2rt}
\alias{ri2rt}
\title{Retention Time Index to Retention Time convertion}
\description{
  Convert retention time indices to retention times indices based on observed FAME RI
  and their standard values.
}
\usage{
ri2rt(riTime, rt.observed, ri.standard)
}
\arguments{
  \item{riTime}{ And RI vector or matrix to convert to Retention Time.}
  \item{rt.observed}{The observed FAME RT's. It could be a vector or a matrix.}
  \item{ri.standard}{ The standard RI for each FAME }
}
\details{
	This function is the inverse of \code{\link{rt2ri}}.
}
\value{
  The converted RT
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{RIcorrect}}, \code{\link{FAMEoutliers}} }
\examples{
  # RI standards
  standard <- c(100, 200, 300, 400, 500)
  # observed standard retention times
  observed <- c(10.4, 19.3, 32.4, 40.2, 50.3)
  # a random set of retention times
  RI       <- runif(100,90,600)	
  # the corrected RIs
  RT       <- ri2rt(RI, observed, standard)
}
