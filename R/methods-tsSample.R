#
# tsSample class methods
#

# valid CDF extensions
._cdf_ext <- c('nc4', 'cdf')

setGeneric("CDFfiles", function(obj) standardGeneric("CDFfiles"))
setMethod("CDFfiles", "tsSample", function(obj) {
    obj@CDFfiles
})

setGeneric("CDFfiles<-", function(obj, value) standardGeneric("CDFfiles<-"))
setReplaceMethod("CDFfiles", "tsSample", function(obj, value) {
    obj@CDFfiles <- value
    validObject(obj)
    obj
})

setGeneric("RIfiles", function(obj) standardGeneric("RIfiles"))
setMethod("RIfiles", "tsSample", function(obj) {
    obj@RIfiles
})

setGeneric("RIfiles<-", function(obj, value) standardGeneric("RIfiles<-"))
setReplaceMethod("RIfiles", "tsSample", function(obj, value) {
    obj@RIfiles <- value
    validObject(obj)
    obj
})

setGeneric("sampleNames", function(obj) standardGeneric("sampleNames"))
setMethod("sampleNames", "tsSample", function(obj) obj@Names )
setGeneric("sampleNames<-", function(obj, value) standardGeneric("sampleNames<-"))
setReplaceMethod("sampleNames", "tsSample", function(obj, value) {
    obj@Names <- value
    # make sure rownames are equal
    rownames(obj@data) <- obj@Names
    validObject(obj)
    obj
})

setGeneric("sampleDays", function(obj) standardGeneric("sampleDays"))
setMethod("sampleDays", "tsSample", function(obj) obj@days )
setGeneric("sampleDays<-", function(obj, value) standardGeneric("sampleDays<-"))
setReplaceMethod("sampleDays", "tsSample", function(obj, value) {
    n <- length(obj@CDFfiles)
    obj@days <- if(length(value) == 1) rep(value, n) else value
    validObject(obj)
    obj
})

setGeneric("CDFpath", function(obj) standardGeneric("CDFpath"))
setMethod("CDFpath", "tsSample", function(obj) .dirname(obj@CDFfiles) )
setGeneric("CDFpath<-", function(obj, value) standardGeneric("CDFpath<-"))
setReplaceMethod("CDFpath", "tsSample", function(obj, value) {
    obj@CDFfiles <- .setpath(obj@CDFfiles, value)
    validObject(obj)
    obj
})

setGeneric("RIpath", function(obj) standardGeneric("RIpath"))
setMethod("RIpath", "tsSample", function(obj) .dirname(obj@RIfiles) )
setGeneric("RIpath<-", function(obj, value) standardGeneric("RIpath<-"))
setReplaceMethod("RIpath", "tsSample", function(obj, value) {
    obj@RIfiles <- .setpath(obj@RIfiles, value)
    validObject(obj)
    obj
})

setGeneric("sampleData", function(obj) standardGeneric("sampleData"))
setMethod("sampleData", "tsSample", function(obj) obj@data )
setGeneric("sampleData<-", function(obj, value) standardGeneric("sampleData<-"))
setReplaceMethod("sampleData", "tsSample", function(obj, value) {
    obj@data <- value
    # make sure rownames are equal
    rownames(obj@data) <- obj@Names
    validObject(obj)
    obj
})

setMethod("length", "tsSample", function(x) length(x@CDFfiles))
setMethod("show", "tsSample", function(object) {
	cat("An object of class 'tsSample':\n")
	cat(" Total samples : ", length(object), "\n")
	cat(" Msrt. days    : ", unique(sampleDays(object)), "\n")
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
	else if(!all(object@Names == rownames(object@data)))
		paste("Invalid object data. Rownames must be equal in metadata")
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
                  x@data <- x@data[i,,drop=FALSE] # must be data.frame
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
          function(.Object, CDFfiles, Names, RIfiles, days, data,
                   CDFpath, RIpath, ftype = c("binary", "text"))
          {
            if(missing(Names))
                Names <- .trim_file_ext(basename(CDFfiles), ._cdf_ext)

            if(missing(RIfiles)) {
                ftype <- match.arg(ftype)
                ext <- c(binary="dat",text="txt")[ftype]
                RIfiles <- .make_RI_files(CDFfiles, ext, ._cdf_ext)
            }
            if(missing(days))
                days <- getDays(basename(CDFfiles))
            else if(length(days) == 1)
                days <- rep(days, length(CDFfiles))

            if(missing(data))
                data <- data.frame(SAMPLE_NAME = Names, stringsAsFactors = FALSE, row.names = Names)
            else
                rownames(data) <- Names

            # update paths
            if(!missing(CDFpath))
                CDFfiles <- .setpath(CDFfiles, CDFpath)
            if(!missing(RIpath))
                RIfiles <- .setpath(RIfiles, RIpath)

            .Object@Names    <- Names
            .Object@CDFfiles <- CDFfiles
            .Object@RIfiles  <- RIfiles
            .Object@CDFpath  <- "" # unused
            .Object@RIpath   <- "" # unused
            .Object@days     <- as.character(days)
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

# warn if tsSample is old
.check_ts_sample <- function(obj, quiet=FALSE)
{
    if(all(c(obj@CDFpath, obj@RIpath) == ""))
        return(TRUE)
    if(!quiet)
        warning("Old `tsSample` object detected. Update it by calling `tsUpdate`")
    return(FALSE)
}

setGeneric("tsUpdate", function(obj) standardGeneric("tsUpdate"))
setMethod("tsUpdate", "tsSample", function(obj) {
    if(.check_ts_sample(obj, TRUE))
        return(obj)
    obj@CDFfiles <- .setpath(obj@CDFfiles, obj@CDFpath)
    obj@RIfiles <- .setpath(obj@RIfiles, obj@RIpath)

    obj@CDFfiles <- obj@CDFfiles
    obj@RIfiles <- obj@RIfiles
    obj@CDFpath <- obj@RIpath <- ""
    rownames(obj@data) <- obj@Names
    validObject(obj)
    obj
})

setGeneric("ncdf4Convert", function(obj, path, ...) standardGeneric("ncdf4Convert"))
setMethod("ncdf4Convert", "tsSample", function(obj, path, ...) {
    if(!.check_ts_sample(obj))
        stop("Object Update required")

    nc4 <- sprintf("%s.nc4", .trim_file_ext(obj@CDFfiles, c('cdf', 'nc4')))
    if(!missing(path))
        nc4 <- .setpath(nc4, path)
    nil <- mapply(ncdf4_convert, obj@CDFfiles, nc4, MoreArgs=list(...))
    obj@CDFfiles <- nc4
    validObject(obj)
    obj
})

# vim: set ts=4 sw=4 et:
