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
\alias{c,tsSample-method}
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
    \item{\code{CDFpath}:}{\code{"character"}, CDF files path. Deprecated. See Notes. }
    \item{\code{RIpath}:}{\code{"character"}, RI file path. Deprecated. See Notes. }
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
    \item{\code{show}}{\code{signature(object = "tsSample")}: the show function. }
    \item{\code{fileFormat}}{\code{signature(obj = "tsSample")}: Sets or gets the RI file format.
        Options are either "binary" or "text". See note below.}
    \item{\code{c}}{Function to combine multiple objects.}
 }
}
\section{Notes}{
    The method \code{fileFormat} only changes the internal information
    of the file type and not the files themselves. To actually
    change the files, use the functions \code{\link{bin2text}}
    and \code{\link{text2bin}}.

    Note that the slot \code{Names} (i.e., the sample names/identifiers) must
    be unique. This allows sample selection by using sample identifiers as well
    as indices. Also, if columns are selected, the output will be either a vector or
    a \code{data.frame} depending on whether one or more columns were selected.
    More over, it is required that the rownames of the \code{data} slot are
    equal to the sample names slots. This is handled internally. If the \code{Names}
    are not provided, then the CDF files are used instead (directories and
    extension are removed). In this case the file names must be unique.

    The slots \code{CDFpath} and \code{RIpath} are deprecated and not used.
    However, the methods to set or get the paths will still work. The file
    paths is stored on the \code{CDFfiles} and \code{RIfiles} slots.

    The \code{c} method requires that the combination of all objects have
    unique identifiers (names), or in other words, the objects cannot share
    identifiers. Duplicated identifiers are not allowed and an error will be thrown.
    This method also requires that the column names of the \code{data} slot
    (see method \code{sampleData}) are the same, because the method calls
    \code{rbind} internally.
}
\author{Alvaro Cuadros-Inostroza, Matthew Hannah, Henning Redestig }
\seealso{
 \code{\link{ImportSamples}}
}
\examples{
showClass("tsSample")

# get a list of CDF files from a directory
require(TargetSearchData)
CDFpath <- tsd_data_path()
cdffiles <- dir(CDFpath, "cdf", full=TRUE)

# define the RI file path
RIpath  <- "."

# create a tsSample object with minimal info
sampleDescription <- new("tsSample", CDFfiles = cdffiles, RIpath = RIpath)

## ## ## ## ##

# Alternatively, the CDF path and CDF file names can be given separately
# (this was the old TargetSearch behavior)
cdffiles <- basename(cdffiles)

# create the sample object
sampleDescription <- new("tsSample", CDFfiles = cdffiles, CDFpath = CDFpath, RIpath = RIpath)

## ## ## ## ##

# More parameters could be defined:
# define the RI files and the RI path
RIfiles <- sub("cdf$", "txt", paste("RI_", cdffiles, sep = ""))
RIpath  <- "."

# get the measurement days (the four first numbers of the cdf files, in this
# example)
days <- substring(cdffiles, 1, 4)

# sample names (must be unique)
smp_names <- paste("Sample", 1:length(sampleDescription), sep = "_")

# add some sample info
smp_data <- data.frame(CDF_FILE =cdffiles, GROUP = gl(5,3))

# create the sample object
sampleDescription <- new("tsSample", Names = smp_names, CDFfiles = cdffiles, CDFpath = CDFpath,
    RIpath = RIpath, days = days, RIfiles = RIfiles, data = smp_data)

# change the file paths (relative to the working path)
CDFpath(sampleDescription) <- "my_cdfs"
RIpath(sampleDescription)  <- "my_RIs"

# change sample Names
sampleNames(sampleDescription) <- sprintf("S\%03d", 1:length(sampleDescription))

## sample subsetting.
# select samples 1, 3 and 5
(sampleset <- sampleDescription[c(1, 3, 5)])

# or use sample IDs
(sampleset <- sampleDescription[c("S001", "S003", "S005")])

## extract columns
# select column 'GROUP'
(group <- sampleDescription$GROUP)
# or
(group <- sampleDescription[, 'GROUP'])

## change the measurement days (sets the same day for all samples)
sampleDays(sampleDescription) <- "1"

## Note: the length of `measurement days` variable must be 1 or equal
## to the number of samples, otherwise an error will be thrown

## objects can be combined with the c() function
samp1 <- sampleDescription[1:5]
samp2 <- sampleDescription[6:10]
samp3 <- sampleDescription[11:15]
samp <- c(samp1, samp2, samp3)
}
\keyword{classes}
