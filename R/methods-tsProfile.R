
# methods for tsMSdata
setGeneric("retIndex", function(obj) standardGeneric("retIndex"))
setMethod("retIndex", "tsMSdata", function(obj) obj@RI)
setGeneric("retIndex<-", function(obj, value) standardGeneric("retIndex<-"))
setReplaceMethod("retIndex", "tsMSdata", function(obj, value) { obj@RI <- value
 validObject(obj)
 obj })

setGeneric("retTime", function(obj) standardGeneric("retTime"))
setMethod("retTime", "tsMSdata", function(obj) obj@RT)
setGeneric("retTime<-", function(obj, value) standardGeneric("retTime<-"))
setReplaceMethod("retTime", "tsMSdata", function(obj, value) { obj@RT <- value
 validObject(obj)
 obj })

setGeneric("Intensity", function(obj) standardGeneric("Intensity"))
setMethod("Intensity", "tsMSdata", function(obj) obj@Intensity)
setGeneric("Intensity<-", function(obj, value) standardGeneric("Intensity<-"))
setReplaceMethod("Intensity", "tsMSdata", function(obj, value) { obj@Intensity <- value
 validObject(obj)
 obj })

setGeneric("as.list.tsMSdata", function(x, ...) standardGeneric("as.list.tsMSdata"))
setMethod("as.list.tsMSdata", "tsMSdata", function(x, ...) {
    am <- function(z) {
        nm <- if(is.null(names(z))) as.character(1:length(z)) else names(z)
        rn <- unlist(sapply(z, rownames))
        cn <- colnames(z[[1]])
        m  <- length(rn)
        n  <- length(cn)
        i  <- rep(1:length(z), sapply(z,nrow))
        y <- matrix(0,m,n, dimnames=list(rn,cn))
        for(k in 1:length(z)) {
            y[k==i,] <- z[[k]]
        }
        attr(y, "index") <- nm[i]
        y
    }
    int <- am(x@Intensity)
    ri  <- am(x@RI)
    if(length(x@RT) > 0) {
        rt  <- am(x@RT)
        return(list(Intensity=int, RI=ri, RT=rt))
    } else {
        return(list(Intensity=int, RI=ri))
    }
})

setMethod("show", "tsMSdata", function(object) {
	cat("An object of class 'tsMSdata':\n")
	cat(" Number of samples: ", ncol(object@RI[[1]]), "\n")
	cat(" Number of metabolites: ", length(object@RI), "\n")
})

validMSdataObject <- function(object) {
    # get lengths of objects
    nri <- length(object@RI)
    nrt <- length(object@RT)
    nin <- length(object@Intensity)

    # get number of samples per component
    nc.ri <- sapply(object@RI, ncol)
    nc.in <- sapply(object@Intensity, ncol)
    nc.rt <- sapply(object@RT, ncol)

    if(nri != nin)
        return(sprintf("Unequal number of `RI` (%d) and `Intensity` (%d) components.", nri, nin))

    if(nri != nrt)
        return(sprintf("Unequal number of `RI` (%d) and `RT` (%d) components.", nri, nrt))

    # check that all the column numbers are the same
    if(any(nc.ri != nc.ri[1]))
        return("Some elements of the `RI` slot have different number of samples.")

    if(any(nc.in != nc.in[1]) || (nc.in[1] != nc.ri[1]))
        return("Some elements of the `Intensity` slot have different number of samples.")

    if(any(nc.rt != nc.rt[1]) || (nc.rt[1] != nc.ri[1]))
        return("Some elements of the `RT` slot have different number of samples.")

    TRUE
}

setValidity("tsMSdata", validMSdataObject)

# methods for tsProfile
setGeneric("profileInfo", function(obj) standardGeneric("profileInfo"))
setMethod("profileInfo", "tsProfile", function(obj) obj@info)
setGeneric("profileInfo<-", function(obj, value) standardGeneric("profileInfo<-"))
setReplaceMethod("profileInfo", "tsProfile", function(obj, value) { obj@info <- value
 validObject(obj)
 obj })

setGeneric("profileInt", function(obj) standardGeneric("profileInt"))
setMethod("profileInt", "tsProfile", function(obj) obj@profInt)
setGeneric("profileInt<-", function(obj, value) standardGeneric("profileInt<-"))
setReplaceMethod("profileInt", "tsProfile", function(obj, value) { obj@profInt <- value
 validObject(obj)
 obj })

setGeneric("profileRI", function(obj) standardGeneric("profileRI"))
setMethod("profileRI", "tsProfile", function(obj) obj@profRI)
setGeneric("profileRI<-", function(obj, value) standardGeneric("profileRI<-"))
setReplaceMethod("profileRI", "tsProfile", function(obj, value) { obj@profRI <- value
 validObject(obj)
 obj })

setGeneric("profileRT", function(obj) standardGeneric("profileRT"))
setMethod("profileRT", "tsProfile", function(obj) obj@profRT)
setGeneric("profileRT<-", function(obj, value) standardGeneric("profileRT<-"))
setReplaceMethod("profileRT", "tsProfile", function(obj, value) { obj@profRT <- value
 validObject(obj)
 obj })

setMethod("show", "tsProfile", function(object) {
	cat("An object of class 'tsMSdata':\n")
	cat(" Profile Information:\n")
	print(head(object@info), 5)
	if(nrow(object@info) > 5) cat("    ", nrow(object@info) - 5,"lines more...\n")
})

setValidity("tsProfile", function(object) {
    validMSdataObject(object)
    n <- length(object@RI)
    s <- sapply(object@RI, ncol)[1]

    if( (nr <- nrow(object@profRI)) != n )
        paste0("Invalid number of rows in slot `profRI`. Found ", nr, ", expected ", n, ".")
    else if( (nc <- ncol(object@profRI)) != s )
        paste0("Invalid number of columns in slot `profRI`. Found ", nc, ", expected ", s, ".")

    else if( (nr <- nrow(object@profRT)) != n )
        paste0("Invalid number of rows in slot `profRT`. Found ", nr, ", expected ", n, ".")
    else if( (nc <- ncol(object@profRT)) != s )
        paste0("Invalid number of columns in slot `profRT`. Found ", nc, ", expected ", s, ".")

    else if( (nr <- nrow(object@profInt)) != n )
        paste0("Invalid number of rows in slot `profInt`. Found ", nr, ", expected ", n, ".")
    else if( (nc <- ncol(object@profInt)) != s )
        paste0("Invalid number of columns in slot `profInt`. Found ", nc, ", expected ", s, ".")

    else if( (nr <- nrow(object@info)) != n )
        paste0("Invalid number of rows in slot `info`. Found ", nr, ", expected ", n, ".")
    else
        TRUE
})

setGeneric("as.list.tsProfile", function(x, ...) standardGeneric("as.list.tsProfile"))
setMethod("as.list.tsProfile", "tsMSdata", function(x, ...) as.list.tsMSdata(x))
