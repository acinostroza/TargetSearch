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
\alias{as.list,tsMSdata-method}
\alias{as.list.tsMSdata,tsMSdata-method}
\alias{as.list.tsMSdata}
\alias{as.list.tsProfile,tsMSdata-method}
\alias{as.list.tsProfile}

\title{Class for representing MS data}
\description{This is a class to represent MS data obtained from the sample.}
\section{Objects from the Class}{
Objects be created by calls of the form 
}
\section{Slots}{
	 \describe{
    \item{\code{RI}:}{\code{"list"}, a list containing an RI matrix, one matrix per metabolite}
    \item{\code{RT}:}{\code{"list"}, a list containing an RT matrix, one matrix per metabolite}
    \item{\code{Intensity}:}{\code{"list"}, a list containing a peak intensity matrix, one matrix per metabolite }
  }
}
\section{Methods}{
  \describe{
    \item{Intensity}{\code{signature(obj = "tsMSdata")}: gets the peak intensity list. }
    \item{Intensity<-}{\code{signature(obj = "tsMSdata")}: gets the peak intensity list. }
    \item{retIndex}{\code{signature(obj = "tsMSdata")}: gets RT list. }
    \item{retIndex<-}{\code{signature(obj = "tsMSdata")}: sets the RI list. }
    \item{retTime}{\code{signature(obj = "tsMSdata")}: gets the RT list. }
    \item{retTime<-}{\code{signature(obj = "tsMSdata")}: sets the RT list. }
    \item{show}{\code{signature(object = "tsMSdata")}: show function. }
    \item{as.list}{\code{signature(object = "tsMSdata")}: coerce a list object. See details }
	 }
}
\details{
  The method \code{as.list} converts every slot (\code{RI}, \code{RT}, and \code{Intensity})
  of a \code{tsMSdata} object into a matrix. The converted matrices
  are stored in a list. Each converted matrix has an attribute called 'index' that
  relates the metabolite index with the respective rows. The components
  of the resulting list are named as the slots. If the slot \code{RT} is
  not defined or empty, then the output list will have only two components.
  ('RT' and 'Intensity').
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }

\seealso{
\code{\link{FindPeaks}}, \code{\link{peakFind}}
}
\examples{
showClass("tsMSdata")
}
\keyword{classes}
