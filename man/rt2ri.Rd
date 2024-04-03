\name{rt2ri}
\alias{rt2ri}
\title{Retention Time to Retention Time Index conversion}
\description{
  Convert retention times to retention indices based on observed FAME RI
  and their standard values.
}
\usage{
rt2ri(rtTime, observed, standard)
}
\arguments{
  \item{rtTime}{ The extracted RT's to convert }
  \item{observed}{The observed FAME RT's }
  \item{standard}{ The standard RI for each FAME }
}
\details{
  Linear interpolation, interpolation outside bounds are done
  with continued linear interpolation from the last two FAME's
}
\value{
  The converted RI
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{RIcorrect}}, \code{\link{FAMEoutliers}} }
\examples{
  # RI standards
  standard <- c(100, 200, 300, 400, 500)
  # observed standard retention times
  observed <- c(10.4, 19.3, 32.4, 40.2, 50.3)
  # a random set of retention times
  RT       <- runif(100,1,60)	
  # the corrected RIs
  RI       <- rt2ri(RT, observed, standard)
}
