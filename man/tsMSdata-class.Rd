\name{tsMSdata-class}
\docType{class}
\alias{tsMSdata-class}
\alias{Intensity,tsMSdata-method}
\alias{Intensity}
\alias{Intensity<-,tsMSdata-method}
\alias{Intensity<-}
\alias{retIndex,tsMSdata-method}
\alias{retIndex}
\alias{retIndex<-,tsMSdata-method}
\alias{retIndex<-}
\alias{retTime,tsMSdata-method}
\alias{retTime}
\alias{retTime<-,tsMSdata-method}
\alias{retTime<-}
\alias{show,tsMSdata-method}

\title{Class for representing MS data}
\description{This is a class to represent MS data obtained from the sample.}
\section{Objects from the Class}{
Objects be created by calls of the form 
}
\section{Slots}{
	 \describe{
    \item{\code{RI}:}{\code{"matrix"}, an RI matrix.}
    \item{\code{RT}:}{\code{"matrix"}, an RT matrix.}
    \item{\code{Intensity}:}{\code{"matrix"}, an peak intensity matrix. }
  }
}
\section{Methods}{
  \describe{
    \item{Intensity}{\code{signature(obj = "tsMSdata")}: gets the peak intensity matrix. }
    \item{Intensity<-}{\code{signature(obj = "tsMSdata")}: gets the peak intensity matrix. }
    \item{retIndex}{\code{signature(obj = "tsMSdata")}: gets RT matrix. }
    \item{retIndex<-}{\code{signature(obj = "tsMSdata")}: sets the RI matrix. }
    \item{retTime}{\code{signature(obj = "tsMSdata")}: gets the RT matrix. }
    \item{retTime<-}{\code{signature(obj = "tsMSdata")}: sets the RT matrix. }
    \item{show}{\code{signature(object = "tsMSdata")}: show function. }
	 }
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }

\seealso{
\code{\link{FindPeaks}}, \code{\link{peakFind}}
}
\examples{
showClass("tsMSdata")
}
\keyword{classes}
