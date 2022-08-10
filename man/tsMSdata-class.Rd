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
  Objects be created by calling the functions \code{\link{FindPeaks}} or
  \code{\link{peakFind}}. See their respective documentation.
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

  These objects are created by the functions \code{\link{FindPeaks}} and
  \code{\link{peakFind}} and there is really no need to create them manually,
  so the examples focus on simple data manipulations.

  The assignment methods \code{retIndex<-}, \code{retTime<-}, and
  \code{Intensity<-} are exposed for convenience, though they should not be
  needed in general. These objects are tightly coupled to the reference library
  (\code{\linkS4class{tsLib}}) and sample (\code{\linkS4class{tsSample}})
  objects, so a modification could cause unexpected errors. Use under your own
  risk if you know what you are doing.
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }

\seealso{
\code{\link{FindPeaks}}, \code{\link{peakFind}}
}
\note{
  For reasons  explained in the details section, the methods \code{retIndex<-},
  \code{retTime<-}, and \code{Intensity<-} are considered \strong{deprecated}
  and could be removed without notice in the feature.
}
\examples{
showClass("tsMSdata")

# some example manipulations using object peakData
data(TSExample)
peakData

# extract the intensities of the metabolites
int <- Intensity(peakData)
# `int` is a list of matrices where rows are m/z values and columns are samples
head(int)

# extract retention indices (same structure for the intensities)
ri <- retIndex(peakData)
head(ri)

# coerce the whole object into a list of matrices (see details)
data <- as.list(peakData)

# it is possible to use the assignment methods `retIndex<-` and `Intensity<-`
# in these objects, but it is not recommended (only if you know what you are
# doing) for example for small tweaks (not sure what this accomplishes however)
ri[[1]][1,1] <- 234000
retIndex(peakData) <- ri
}
\keyword{classes}
