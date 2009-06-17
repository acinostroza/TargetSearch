\name{ImportFameSettings}
\alias{ImportFameSettings}
\title{ Retention time markers settings }
\description{
  This function imports a list of retention standard markers.
}
\usage{
ImportFameSettings(tmp.file = NA, mass = NA, ...)
}
\arguments{
  \item{tmp.file}{ A character string naming a file with standard markers. }
  \item{mass}{ The m/z standard markers. }
  \item{\dots}{ Other options passed to \code{read.delim} function. }
}
\details{
    The standard marker file is a tab-delimited text file with 3 or 4 columns.
    Column names doesn't matter, however they must be in the following order.
 
  \itemize{
  \item \code{LowerLimit} - The Retention time lower limit in seconds.
  \item \code{UpperLimit} - The Retention time upper limit in seconds.
  \item \code{RIstandard} - The RI value of that standard.
  \item \code{mzMarker}   - The m/z standard markers. This is optional
  since it could be specified with the \code{mass} parameter.
  }
    
  If \code{mass} is a numeric vector of length equal to 1, this value will
  be applied to all markers. Otherwise, the length of the vector must be
  equal to the number of standards.
    If no arguments are given, a default object will be returned.
}
\value{
    A \code{tsRim} object.
}
\examples{
# get the RI marker definition file
cdfpath  <- file.path(.find.package("TargetSearchData"), "gc-ms-data")
rim.file <- file.path(cdfpath, "rimLimits.txt")

# set the mass marker to 87
mass <- 87

# load the definition
rimLimits <- ImportFameSettings(rim.file, mass = mass)

# sometimes you need to change the limits of a particular standard
rimLimits(rimLimits)[2,] <- c(410, 450)

# to change the mass value
rimMass(rimLimits) <- 85

}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{ \code{\link{RIcorrect}}, \code{\linkS4class{tsRim}} }
