\name{tsRim-class}
\docType{class}
\alias{tsRim-class}
\alias{[,tsRim-method}
\alias{rimLimits,tsRim-method}
\alias{rimLimits}
\alias{rimLimits<-,tsRim-method}
\alias{rimLimits<-}
\alias{rimMass,tsRim-method}
\alias{rimMass}
\alias{rimMass<-,tsRim-method}
\alias{rimMass<-}
\alias{rimStandard,tsRim-method}
\alias{rimStandard}
\alias{rimStandard<-,tsRim-method}
\alias{rimStandard<-}
\alias{length,tsRim-method}
\alias{c,tsRim-method}

\title{Class for representing retention index markers}
\description{This is a class to represent retention index markers.}
\section{Objects from the Class}{
Objects can be created by the function \code{\link{ImportFameSettings}} or by calls
of the form \code{new("tsRim", limits = [two column matrix with time limits],
standard = [a vector with RI standards], mass = [m/z marker])}.

}
\section{Slots}{
  \describe{
    \item{\code{limits}:}{\code{"matrix"}, two column matrix with lower and upper
        limits where the standards will be search. One row per standard. }
    \item{\code{standard}:}{\code{"numeric"}, the marker RI values. }
    \item{\code{mass}:}{\code{"numeric"}, the m/z marker. }
  }
}
\section{Methods}{
  \describe{
    \item{\code{[}}{\code{signature(x = "tsRim")}: Selects a subset of markers.}
    \item{\code{rimLimits}}{\code{signature(obj = "tsRim")}: gets the time limits. }
    \item{\code{rimLimits<-}}{\code{signature(obj = "tsRim")}: sets the time limits. }
    \item{\code{rimMass}}{\code{signature(obj = "tsRim")}: gets the m/z marker. }
    \item{\code{rimMass<-}}{\code{signature(obj = "tsRim")}: sets the m/z marker. }
    \item{\code{rimStandard}}{\code{signature(obj = "tsRim")}: gets the standars. }
    \item{\code{rimStandard<-}}{\code{signature(obj = "tsRim")}: sets the standars. }
    \item{\code{length}}{Returns the number of retention index standards.}
    \item{\code{c}}{Function to combine multiple objects.}
  }
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{
 \code{\link{ImportFameSettings}}
}
\examples{
showClass("tsRim")

# create a rimLimit object:
#  - set the lower (first column) and upper (second column) time limites to
#    search for standards.
Lim <- rbind(c(200, 300), c(400,450), c(600,650))
#  - set the retention indices of the standard
Std <- c(250000, 420000, 630000)
#  - set the mass marker
mass <- 87
#  - create the object
rimLimits <- new("tsRim", limits = Lim, standard = Std, mass = mass)
#  - get the number of standards using `length`
length(rimLimits)

# sometimes you need to change the limits of a particular standard
rimLimits(rimLimits)[2,] <- c(410, 450)

# to change the mass value
rimMass(rimLimits) <- 85

# to select a subset
rim <- rimLimits[1:2]

# remove a marker (number 3 in this case)
rim <- rimLimits[-3]

# objects can be combined using the `c` method
obj <- new("tsRim", limits = cbind(100 * 1:6, 100 * 1:6 + 20),
           standard = 1000 * 1:6, mass = 100) # a random object
# combine objects
obj <- c(rimLimits, obj)
}
\keyword{classes}
% vim: set ts=4 sw=4 et:
