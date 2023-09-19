\name{tsProfile-class}
\Rdversion{1.1}
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
\description{
  This is class used to represent a metabolite profile obtained with
  \code{\link{TargetSearch}}. This class contains the data you are most
  likely to use in downstream data analysis.
}
\section{Objects from the Class}{
  Objects of this class are created by the functions \code{\link{Profile}} or
  \code{\link{ProfileCleanUp}}. There is no need to generate objects using
  \code{\link[methods]{new}} methods.
}
\section{Slots}{
  \describe{
    \item{\code{info}:}{\code{"data.frame"}, the metabolite profile information.
      See details.}
    \item{\code{RI}:}{\code{"list"}, a list of RI matrices, one matrix per
      metabolite, in which rows represent select/top masses and columns
      correspond to samples.}
    \item{\code{RT}:}{\code{"list"}, a list of RT matrices, one matrix per
      metabolite, in which rows represent select/top masses and columns
      correspond to samples. }
    \item{\code{Intensity}:}{\code{"list"}, a list of peak-intensity matrices,
      one matrix per metabolite, and within each matrix, rows represent
      select/top masses and columns correspond to samples.}
    \item{\code{profRI}:}{\code{"matrix"}, the profile RI matrix.}
    \item{\code{profRT}:}{\code{"matrix"}, the profile RT matrix.}
    \item{\code{profInt}:}{\code{"matrix"}, the profile Intensity matrix.}
  }
}
\section{Extends}{
  This class extends \code{\linkS4class{tsMSdata}} directly. Methods that
  apply to that class apply to \code{\linkS4class{tsProfile}}.
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
\details{
  The metabolite profile information (slot \code{info}) is a \code{data.frame}
  that contains information regarding each target metabolite in the reference
  library (object \code{\linkS4class{tsLib}}). The information can be accessed
  with the method \code{profileInfo}. The profile include the following
  columns,

  \itemize{
    \item \code{Name} The metabolite/analyte name.
    \item \code{Lib_RI} The reference or expected RI (retention index).
    \item \code{Mass_count} The number of correlating masses.
    \item \code{Non_consecutive_Mass_count} Same as mass count, but
      discounting consecutive masses that correspond with isotopic peaks, which
      often correlate anyway.
    \item \code{Sample_Count_per_Mass} The number of samples in which each
      correlating mass was found. The values are separated by colons (\code{;}).
      The values of the actual masses are in following column.
    \item \code{Masses} The correlating masses separated by colons (\code{;}).
    \item \code{RI} The average RI detected in the samples.
    \item \code{Score_all_masses} The similarity score calculated using the
      average intensity of all the masses that were searched for, regardless of
      whether they are correlating masses.
    \item \code{Score_cor_masses} Same as above, but only correlating masses
      are considered.
  }

  Other columns could be present by they are carried over from the metabolite
  library \code{\linkS4class{tsLib}} object.

  If the object is created by the function \code{\link{ProfileCleanUp}} instead,
  additional columns will be created with the prefix \code{Cor_} followed
  by a string, such as \code{Cor_RI}. These columns contain information
  of \dQuote{collapsed} metabolites that correlate to mean that they could
  conflict. The collapsed values are often separated by a vertical bar \code{|}.

  The slots \code{profRT}, \code{profRI} and \code{profInt} contain averaged
  (median) values for each the matrices in the slots \code{RT}, \code{RI},
  and \code{Intensity}. The values are computed only considering the masses
  that correlate (the other masses are ignored). They represent an average of
  RT, RI and Intensity per metabolite per sample, given that each selected mass
  has their own values and often only and average is needed. These matrices are
  accessed with the methods \code{profileRT}, \code{profileRI} and
  \code{profileRT}.

  Similar to the \code{\linkS4class{tsMSdata}} assignment methods, the
  \dQuote{setters} methods are exposed for convenience, though they should
  not be required. The reason is that all the objects are tightly tied and
  unexpected modifications could result in unintended errors. Use only if you
  know what you are doing.
}
\author{Alvaro Cuadros-Inostroza}
\seealso{
  \code{\link{Profile}}, \code{\link{ProfileCleanUp}}, \code{\linkS4class{tsMSdata}}
}
\note{
  For reasons explained in the details section, the methods \code{profileRI<-},
  \code{profileRT<-}, \code{profileInt<-}, and \code{profileInfo} are considered
  \strong{deprecated} and could be removed without notice in the feature.
}
\examples{
showClass("tsProfile")

# some example manipulations using object peakData
data(TSExample)
metabProfile

# extract the profile information of the metabolite analysis
info <- profileInfo(metabProfile)
# `info` is a data.frame with the information described in details
head(info)

# extract profiled metabolite intensities
mat <- profileInt(metabProfile)
head(mat)

# It is possible to use the assignment methods such as `retIndex<-`,
# `Intensity<-`, `profileInt`, etc., in these objects, but it isn't
# recommended. Use only if if you know what you are doing.

\dontrun{
  profileInt(metabProfile) <- updatedMetabMatrix
  profileInfo(metabProfile) <- updatedMetabInfo
}

# in the mock examples above, the objects updated* are updated data.frame or
# matrices with metabolite results.

}
\keyword{classes}
