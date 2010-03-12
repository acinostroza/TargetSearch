\name{tsProfile-class}
\docType{class}
\alias{tsProfile-class}
\alias{profileInfo,tsProfile-method}
\alias{profileInfo}
\alias{profileInfo<-,tsProfile-method}
\alias{profileInfo<-}
\alias{profileInt,tsProfile-method}
\alias{profileInt}
\alias{profileInt<-,tsProfile-method}
\alias{profileInt<-}
\alias{profileRI,tsProfile-method}
\alias{profileRI}
\alias{profileRI<-,tsProfile-method}
\alias{profileRI<-}
\alias{profileRT,tsProfile-method}
\alias{profileRT}
\alias{profileRT<-,tsProfile-method}
\alias{profileRT<-}
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
    \item{\code{RI}:}{\code{"list"}, a list containing an RI matrix, one matrix per metabolite}
    \item{\code{RT}:}{\code{"list"}, a list containing an RT matrix, one matrix per metabolite}
    \item{\code{Intensity}:}{\code{"list"}, a list containing a peak intensity matrix, one matrix per metabolite }
    \item{\code{profRI}:}{\code{"matrix"}, the profile RI matrix.}
    \item{\code{profRT}:}{\code{"matrix"}, the profile RT matrix.}
    \item{\code{profInt}:}{\code{"matrix"}, the profile Intensity matrix.}
  }
}
\section{Extends}{
Class \code{\linkS4class{tsMSdata}}, directly.
}
\section{Methods}{
  \describe{
    \item{profileInfo}{\code{signature(obj = "tsProfile")}: get the profile information. }
    \item{profileInfo<-}{\code{signature(obj = "tsProfile")}: set the profile information. }
    \item{profileInt}{\code{signature(obj = "tsProfile")}: get the profile intensity matrix. }
    \item{profileInt<-}{\code{signature(obj = "tsProfile")}: set the profile intensity matrix. }
    \item{profileRI}{\code{signature(obj = "tsProfile")}: get the profile RI matrix. }
    \item{profileRI<-}{\code{signature(obj = "tsProfile")}: set the profile RI matrix. }
    \item{profileRT}{\code{signature(obj = "tsProfile")}: get the profile RT matrix. }
    \item{profileRT<-}{\code{signature(obj = "tsProfile")}: set the profile RT matrix. }
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
