\name{makeIndex}
\alias{makeIndex}
\title{Make an index vector out of selective/top masses}
\description{
  A function that makes an index for selective or top masses.
  Replaces the method \code{\link{libId}}.
}
\usage{
  makeIndex(lib, sel=TRUE)
}
\arguments{
  \item{lib}{A \linkS4class{tsLib} object.}
  \item{sel}{A logical flag. If \code{TRUE} return an integer index for
    selective masses, otherwise for top masses.}
}
\value{
  An integer vector, where each number is associated with an element of
  the library object, and its multiplicity with the number of selective
  or top masses.
}
\details{
  This is an internal function that should not be used. It is exposed
  because the method it replaces, \code{\link{libId}}, was originally
  exposed, so compatibility is not lost.
}
\note{
  It is equivalent to the deprecated method \code{\link{libId}}.
}
\examples{
# load example objects
data(TSExample)
# generate an index for selective masses
(idx <- makeIndex(refLibrary, TRUE))
}
\author{Alvaro Cuadros-Inostroza}
\seealso{ \code{\linkS4class{tsLib}} }

% vim: set ts=2 sw=2 et:
