\name{tsLib-class}
\docType{class}
\alias{tsLib-class}
\alias{[,tsLib-method}
\alias{$,tsLib-method}
\alias{libId,tsLib-method}
\alias{libId}
\alias{length,tsLib-method}
\alias{libData,tsLib-method}
\alias{libData}
\alias{libName,tsLib-method}
\alias{libName}
\alias{libName<-,tsLib-method}
\alias{libName<-}
\alias{libRI,tsLib-method}
\alias{libRI}
\alias{libRI<-,tsLib-method}
\alias{libRI<-}
\alias{medRI,tsLib-method}
\alias{medRI}
\alias{medRI<-,tsLib-method}
\alias{medRI<-}
\alias{refLib,tsLib-method}
\alias{refLib}
\alias{RIdev,tsLib-method}
\alias{RIdev}
\alias{RIdev<-,tsLib-method}
\alias{RIdev<-}
\alias{quantMass,tsLib-method}
\alias{quantMass}
\alias{quantMass<-,tsLib-method}
\alias{quantMass<-}
\alias{selMass,tsLib-method}
\alias{selMass}
\alias{selMass<-,tsLib-method}
\alias{selMass<-}
\alias{show,tsLib-method}
\alias{spectra,tsLib-method}
\alias{spectra}
\alias{spectra<-,tsLib-method}
\alias{spectra<-}
\alias{topMass,tsLib-method}
\alias{topMass}
\alias{topMass<-,tsLib-method}
\alias{topMass<-}

\title{Class for representing a reference library }
\description{ This is a class representation of a reference library.	}
\section{Objects from the Class}{
Objects can be created by the function \code{\link{ImportLibrary}}.
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{\code{"character"}, the metabolite or analyte names.}
    \item{\code{RI}:}{\code{"numeric"}, the expected retention time indices (RI) of the metabolites/analytes.}
    \item{\code{medRI}:}{\code{"numeric"}, the median RI calculated from the samples.}
    \item{\code{RIdev}:}{\code{"matrix"}, the RI deviation windows, k = 1,2,3. A three column matrix}
    \item{\code{selMass}:}{\code{"list"}, every component is a numeric vector containing the selective masses. }
    \item{\code{topMass}:}{\code{"list"}, every component is a numeric vector containing the top masses. }
    \item{\code{quantMass}:}{\code{"numeric"}, the mass used for quantification.}
    \item{\code{libData}:}{\code{"data.frame"}, additional library information. }
    \item{\code{spectra}:}{\code{"list"}, the metabolite spectra. Each component is a two column matrix: m/z and intensity. }
  }
}
\section{Methods}{
  \describe{
    \item{\code{[}}{\code{signature(x = "tsLib")}: Selects a subset of metabolites from the library.}
    \item{\code{$name}}{\code{signature(x = "tsLib")}: Access column \code{name} of \code{libData} slot. }
    \item{\code{libId}}{\code{signature(obj = "tsLib")}: Returns a vector of indices. }
    \item{\code{length}}{\code{signature(x = "tsLib")}: returns the length of the library. i.e., number of metabolites.}
    \item{\code{libData}}{\code{signature(obj = "tsLib")}: gets the \code{libData} slot.}
    \item{\code{libName}}{\code{signature(obj = "tsLib")}: gets the \code{Name} slot. }
    \item{\code{libRI}}{\code{signature(obj = "tsLib")}: gets the \code{RI} slot.  }
    \item{\code{medRI}}{\code{signature(obj = "tsLib")}: gets the \code{medRI} slot.  }
    \item{\code{refLib}}{\code{signature(obj = "tsLib")}: Low level method to create a matrix representation of the library.}
    \item{\code{RIdev}}{\code{signature(obj = "tsLib")}: gets the RI deviations. }
    \item{\code{RIdev<-}}{\code{signature(obj = "tsLib")}: sets the RI deviations. }
    \item{\code{quantMass}}{\code{signature(obj = "tsLib")}: gets the quantification mass. }
    \item{\code{quantMass<-}}{\code{signature(obj = "tsLib")}: sets the quantification mass. }
    \item{\code{selMass}}{\code{signature(obj = "tsLib")}: gets the selective masses. }
    \item{\code{show}}{\code{signature(object = "tsLib")}: show method. }
    \item{\code{spectra}}{\code{signature(obj = "tsLib")}: gets the spectra. }
    \item{\code{topMass}}{\code{signature(obj = "tsLib")}: gets the top masses. }
	 }
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{
	\code{\link{ImportLibrary}}
}
\examples{
showClass("tsLib")

# define some metabolite names
libNames   <- c("Metab1", "Metab2", "Metab3")
# the expected retention index
RI         <- c(100,200,300)
# selective masses to search for. A list of vectors.
selMasses  <- list(c(95,204,361), c(87,116,190), c(158,201,219))
# define the retention time windows to look for the given selective masses.
RIdev      <- matrix(rep(c(10,5,2), length(libNames)), ncol = 3, byrow = TRUE)
# Set the mass spectra. A list object of two-column matrices, or set to 
# an empty list if the spectra is not available
spectra    <- list()
# some extra information about the library
libData    <- data.frame(Name = libNames, Lib_RI = RI)
# create a reference library object
refLibrary <- new("tsLib", Name = libNames, RI = RI, medRI = RI, RIdev = RIdev,
			selMass = selMasses, topMass = selMasses, spectra = spectra, libData = libData) 

# get the metabolite names
libName(refLibrary)
# set new names
libName(refLibrary) <- c("Metab01", "Metab02", "Metab03")

# get the expected retention times
libRI(refLibrary)
# set the retention time index for metabolite 3 to 310 seconds
libRI(refLibrary)[3] <- 310
# change the seleccion and top masses of metabolite 3
selMass(refLibrary)[[3]] <- c(158,201,219,220,323)
topMass(refLibrary)[[3]] <- c(158,201,219,220,323)
# change the retention time deviations
RIdev(refLibrary)[3,] <- c(8,4,1)

}
\keyword{classes}
