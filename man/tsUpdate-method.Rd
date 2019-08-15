\docType{methods}
\name{tsUpdate,tsSample-method}
\alias{tsUpdate,tsSample-method}
\alias{tsUpdate}
\title{ Methods for Updating TargetSearch objects \code{tsUpdate} }
\description{
    \code{tsUpdate} is a generic function which can be used to update
    and old \sQuote{TargetSearch} class definition. Currently, this
    function is only implemented for \code{\linkS4class{tsSample}} objects.
}
\section{Methods}{
\describe{
    \item{\code{signature(obj = "tsSample")}}{
        Method to update an old \code{\linkS4class{tsSample}} object.
	A change was introduced starting from \sQuote{TargetSearch} version 1.42.0.
}
}}
\examples{
\dontrun{
newObject <- tsUpdate(oldObject)
}
}
\keyword{methods}
\keyword{tsSample}
\author{Alvaro Cuadros-Inostroza}
