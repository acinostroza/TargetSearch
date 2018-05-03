#
# tsSample class methods
#
setGeneric("CDFfiles", function(obj) standardGeneric("CDFfiles"))
setMethod("CDFfiles", "tsSample", function(obj) {
    if(obj@CDFpath != ".")
        file.path(obj@CDFpath,obj@CDFfiles)
    else
        obj@CDFfiles
  })

setGeneric("CDFfiles<-", function(obj, value) standardGeneric("CDFfiles<-"))
setReplaceMethod("CDFfiles", "tsSample", function(obj, value) {
	obj@CDFfiles <- basename(value)
	validObject(obj)
	obj
})

setGeneric("RIfiles", function(obj) standardGeneric("RIfiles"))
setMethod("RIfiles", "tsSample", function(obj) {
	if(obj@RIpath != ".")
		file.path(obj@RIpath,obj@RIfiles)
	else
	obj@RIfiles
})

setGeneric("RIfiles<-", function(obj, value) standardGeneric("RIfiles<-"))
setReplaceMethod("RIfiles", "tsSample", function(obj, value) {
	obj@RIfiles <- basename(value)
	validObject(obj)
	obj
})

setGeneric("sampleNames", function(obj) standardGeneric("sampleNames"))
setMethod("sampleNames", "tsSample", function(obj) obj@Names )
setGeneric("sampleNames<-", function(obj, value) standardGeneric("sampleNames<-"))
setReplaceMethod("sampleNames", "tsSample", function(obj, value) {
	obj@Names <- value
	validObject(obj)
	obj
})

setGeneric("sampleDays", function(obj) standardGeneric("sampleDays"))
setMethod("sampleDays", "tsSample", function(obj) obj@days )
setGeneric("sampleDays<-", function(obj, value) standardGeneric("sampleDays<-"))
setReplaceMethod("sampleDays", "tsSample", function(obj, value) {
	obj@days <- value
	validObject(obj)
	obj
})

setGeneric("CDFpath", function(obj) standardGeneric("CDFpath"))
setMethod("CDFpath", "tsSample", function(obj) obj@CDFpath )
setGeneric("CDFpath<-", function(obj, value) standardGeneric("CDFpath<-"))
setReplaceMethod("CDFpath", "tsSample", function(obj, value) {
	obj@CDFpath <- value
	obj
})

setGeneric("RIpath", function(obj) standardGeneric("RIpath"))
setMethod("RIpath", "tsSample", function(obj) obj@RIpath )
setGeneric("RIpath<-", function(obj, value) standardGeneric("RIpath<-"))
setReplaceMethod("RIpath", "tsSample", function(obj, value) {
	obj@RIpath <- value
	obj
})

setGeneric("sampleData", function(obj) standardGeneric("sampleData"))
setMethod("sampleData", "tsSample", function(obj) obj@data )
setGeneric("sampleData<-", function(obj, value) standardGeneric("sampleData<-"))
setReplaceMethod("sampleData", "tsSample", function(obj, value) {
	obj@data <- value
	validObject(obj)
	obj
})


setMethod("length", "tsSample", function(x) length(x@CDFfiles))
setMethod("show", "tsSample", function(object) {
	cat("An object of class 'tsSample':\n")
	cat(" Number of samples:  ", length(object), "\n")
	cat(" CDF files directory:", object@CDFpath, "\n")
	cat(" RI files directory: ", object@RIpath, "\n")
	cat(" Measurement days:   ", unique(sampleDays(object)), "\n")
	cat("\nSample Data:\n")
	print(head(sampleData(object), 5))
	if(length(object) > 5) cat("    ", length(object) - 5,"lines more...\n")
})

setValidity("tsSample", function(object) {
	n <- length(object@CDFfiles)
	if(length(object@RIfiles) != n)
		paste("Unequal number of CDF and RI files: ", n,", ", length(object@RIfiles), sep = "")
	else if(length(object@Names) != n)
		paste("Unequal number of CDF files and Names: ", n,", ", length(object@Names), sep = "")
	else if(length(object@days) != n)
		paste("Unequal number of CDF files and Days: ", n,", ", length(object@days), sep = "")
	else if(nrow(object@data) != n)
		paste("Unequal number of CDF files and Sample Data: ", n,", ", nrow(object@data), sep = "")
	else if(any(duplicated(object@Names)))
		paste("sample names must be unique. duplicated names found")
	else TRUE
})

setMethod("[", "tsSample",
          function(x, i, j, ..., drop=TRUE)
          {
              if(missing(i) & missing(j)) {
                  validObject(x)
                  return(x)
              }

              if(missing(j)) {
                  if(is.character(i))
                      i <- match(i, x@Names)
                  x@Names <- x@Names[i]
                  x@CDFfiles <- x@CDFfiles[i]
                  x@RIfiles <- x@RIfiles[i]
                  x@days <- x@days[i]
                  x@data <- x@data[i,]
                  validObject(x)
                  return(x)
              } else {
                  x@data[i, j, drop=drop]
              }
          })

setMethod("$", "tsSample", function(x, name) {
    ret <- eval(substitute(sampleData(x)$NAME_ARG, list(NAME_ARG=name)))
    if(is.null(ret))
        warning("Column `", name, "` does not exist")
    ret
})

setMethod("initialize",
          "tsSample",
          function(.Object, Names    = character(0), CDFfiles = character(0),
                            RIfiles  = character(0), CDFpath  = ".", RIpath   = ".",
                            days     = character(0), data     = data.frame(),
                            ftype    = "binary") {
            if (length(CDFfiles) > 0) {
                if(length(Names) == 0)
                    Names <- gsub(".cdf", "", CDFfiles, ignore.case = T)
                if(length(RIfiles) == 0) {
                    ftype <- pmatch(ftype, c("binary","text"))[1]
                    stopifnot(!is.na(ftype))
                    ext <- c("dat","txt")[ftype]
                    RIfiles <- sub("cdf$", ext, paste("RI_", CDFfiles, sep = ""), ignore.case = T)
                }
                if(length(days) == 0)
                    days <- getDays(CDFfiles)
                if(all(dim(data) == 0))
                    data <- data.frame(Names = Names, CDF_FILE = CDFfiles,
                                       MEASUREMENT_DAY = days, RI_FILE = RIfiles,
                                       stringsAsFactors = FALSE)
            }

            .Object@Names    <- Names
            .Object@CDFfiles <- CDFfiles
            .Object@RIfiles  <- RIfiles
            .Object@CDFpath  <- CDFpath
            .Object@RIpath   <- RIpath
            .Object@days     <- days
            .Object@data     <- data
            validObject(.Object)
            .Object
          })

# new methods for 'binary' or 'text' RI-file format.
# the file type is determined by the file extension.
#   *.dat => binary, *.txt => text
# returns either 'binary' or 'text'
setGeneric("fileFormat", function(obj) standardGeneric("fileFormat"))
setMethod("fileFormat", "tsSample", function(obj) {
	tmp <- unique(substring(tolower(obj@RIfiles), nchar(obj@RIfiles)-2))
	if(length(tmp)!=1)
		warning("Incorrect RIfile extensions. Expected '*.dat' or '*.txt'")
	c("binary","text")[pmatch(tmp[1], c("dat", "txt"))]
})

setGeneric("fileFormat<-", function(obj, value) standardGeneric("fileFormat<-"))
setReplaceMethod("fileFormat", "tsSample", function(obj, value) {
	value <- pmatch(value[1], c("binary","text"))
	if(is.na(value))
		stop("Incorrect file format: It should be either 'binary' or 'text'")
	obj@RIfiles <- sub("\\.\\w{3}$", c(".dat",".txt")[value], obj@RIfiles, perl=TRUE)
	if(!is.null(obj@data$RI_FILE))
		obj@data$RI_FILE <- obj@RIfiles
	validObject(obj)
	obj
})

# vim: set ts=4 sw=4 et:
