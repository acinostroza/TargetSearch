\name{tsSample-class}
\docType{class}
\alias{tsSample-class}
\alias{[,tsSample-method}
\alias{$,tsSample-method}
\alias{CDFfiles,tsSample-method}
\alias{CDFfiles}
\alias{CDFfiles<-,tsSample-method}
\alias{CDFfiles<-}
\alias{RIfiles,tsSample-method}
\alias{RIfiles}
\alias{RIfiles<-,tsSample-method}
\alias{RIfiles<-}
\alias{CDFpath,tsSample-method}
\alias{CDFpath}
\alias{CDFpath<-,tsSample-method}
\alias{CDFpath<-}
\alias{RIpath,tsSample-method}
\alias{RIpath}
\alias{RIpath<-,tsSample-method}
\alias{RIpath<-}
\alias{length,tsSample-method}
\alias{sampleData,tsSample-method}
\alias{sampleData}
\alias{sampleData<-,tsSample-method}
\alias{sampleData<-}
\alias{sampleDays,tsSample-method}
\alias{sampleDays}
\alias{sampleDays<-,tsSample-method}
\alias{sampleDays<-}
\alias{sampleNames,tsSample-method}
\alias{sampleNames}
\alias{sampleNames<-,tsSample-method}
\alias{sampleNames<-}
\alias{fileFormat,tsSample-method}
\alias{fileFormat}
\alias{fileFormat<-,tsSample-method}
\alias{fileFormat<-}
\alias{show,tsSample-method}

\title{Class for representing samples}
\description{This is a class to represent a set of samples.}
\section{Objects from the Class}{
Objects can be created by the function \code{\link{ImportSamples}} or by calling
the object generator function.

\code{
new("tsSample", Names = [sample names], CDFfiles = [list of CDF file names],
     RIfiles = [list of RI file names], CDFpath = [CDF files path],
     RIpath  = [RI files path], days = [measurement days],
     data    = [additional sample information], ftype = [RI file format])
}}

\section{Slots}{
  \describe{
    \item{\code{Names}:}{\code{"character"}, the sample names.}
    \item{\code{CDFfiles}:}{\code{"character"}, the list of CDF file names. }
    \item{\code{RIfiles}:}{\code{"character"}, the list of RI file names. }
    \item{\code{CDFpath}:}{\code{"character"}, CDF files path. }
    \item{\code{RIpath}:}{\code{"character"}, RI file path. }
    \item{\code{days}:}{\code{"character"}, measurement days.}
    \item{\code{data}:}{\code{"data.frame"}, additional sample information.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{[}}{\code{signature(x = "tsSample")}: Selects a subset of samples.}
    \item{\code{$name}}{\code{signature(x = "tsSample")}: Access column \code{name} of \code{sampleData} slot.}
    \item{\code{CDFfiles}}{\code{signature(obj = "tsSample")}: list of CDF files. }
    \item{\code{RIfiles}}{\code{signature(obj = "tsSample")}: list of RI files. }
    \item{\code{RIpath}}{\code{signature(obj = "tsSample")}: The RI file path. }
    \item{\code{CDFpath}}{\code{signature(obj = "tsSample")}: The CDF file path. }
    \item{\code{length}}{\code{signature(x = "tsSample")}: number of samples.}
    \item{\code{sampleData}}{\code{signature(obj = "tsSample")}: additional sample information. }
    \item{\code{sampleDays}}{\code{signature(obj = "tsSample")}: measurement days. }
    \item{\code{sampleNames}}{\code{signature(obj = "tsSample")}: sample names. The names must be unique}
    \item{\code{show}}{\code{signature(object = "tsSample")}: the show funtion. }
    \item{\code{fileFormat}}{\code{signature}{obj = "tsSample"}: Sets or gets the RI file format.
        Options are either "binary" or "text". See note below.}
 }
}
\section{Notes}{
    The method \code{fileFormat} only changes the internal information
    of the file type and not the files themselves. To actually
    change the files, use the functions \code{\link{bin2text}}
    and \code{\link{text2bin}}.
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{
 \code{\link{ImportSamples}}
}
\examples{
showClass("tsSample")

# get a list of CDF files from a directory
require(TargetSearchData)
CDFpath <- system.file("gc-ms-data", package = "TargetSearchData")
cdffiles <- dir(CDFpath, "cdf")

# define the RI file path 
RIpath  <- "."

# create the sample object
sampleDescription <- new("tsSample", CDFfiles = cdffiles, CDFpath = CDFpath, RIpath = RIpath) 

##

# More parameters could be defined:
# define the RI files and the RI path
RIfiles <- sub("cdf$", "txt", paste("RI_", cdffiles, sep = ""))
RIpath  <- "."

# get the measurement days (the four first numbers of the cdf files, in this
# example)
days <- substring(cdffiles, 1, 4)

# sample names 
smp_names <- sub(".cdf", "", cdffiles, fixed = TRUE)

# add some sample info
smp_data <- data.frame(CDF_FILE =cdffiles, GROUP = gl(5,3))

# create the sample object
sampleDescription <- new("tsSample", Names = smp_names, CDFfiles = cdffiles, CDFpath = CDFpath,
    RIpath = RIpath, days = days, RIfiles = RIfiles, data = smp_data)

# changing the sample names
sampleNames(sampleDescription) <- paste("Sample", 1:length(sampleDescription), sep = "_")

# changing the file paths (relative to the working path)
CDFpath(sampleDescription) <- "my_cdfs/"
RIpath(sampleDescription)  <- "my_RIs/"

}
\keyword{classes}
