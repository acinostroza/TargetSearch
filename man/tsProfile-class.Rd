\name{tsProfile-class}
\docType{class}
\alias{tsProfile-class}
\alias{profileInfo,tsProfile-method}
\alias{profileInfo}
\alias{profileInfo<-,tsProfile-method}
\alias{profileInfo<-}
\alias{show,tsProfile-method}

\title{Class for representing a MS profile}
\description{This class is to represent a MS profile}
\section{Objects from the Class}{
Objects can be created by the function \code{\link{Profile}} or by

\code{
new("tsMSdata", RI = [retention time index matrix], RT = [retention time matrix],
                Intensity = [peak intensity])
}

}
\section{Slots}{
	 \describe{
    \item{\code{info}:}{\code{"data.frame"}, the profile information. }
    \item{\code{RI}:}{\code{"matrix"}, an RI matrix.}
    \item{\code{RT}:}{\code{"matrix"}, an RT matrix.}
    \item{\code{Intensity}:}{\code{"matrix"}, an peak intensity matrix. }
  }
}
\section{Extends}{
Class \code{\linkS4class{tsMSdata}}, directly.
}
\section{Methods}{
  \describe{
    \item{profileInfo}{\code{signature(obj = "tsProfile")}: get the profile information. }
    \item{profileInfo<-}{\code{signature(obj = "tsProfile")}: set the profile information. }
    \item{show}{\code{signature(object = "tsProfile")}: the show function. }
	 }
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{
	\code{\link{Profile}}, 	\code{\link{ProfileCleanUp}}, 
	\code{\linkS4class{tsMSdata}}
}
\examples{
showClass("tsProfile")
}
\keyword{classes}
