
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
    nri <- length(object@RI)
    nrt <- length(object@RT)
    nin <- length(object@Intensity)
    if(nri != nin)
        return(paste("Unequal number of RIs and Intensities: ", nri,", ", nin, sep = ""))

    # check that all the column numbers are the same
    nc.ri <- sapply(object@RI, ncol)
    if(!all(nc.ri == nc.ri[1]))
        return(paste("Some elements of the RI slot have different number of samples"))
    nc.in <- sapply(object@Intensity, ncol)
    if(!all(nc.in == nc.in[1]))
        return(paste("Some elements of the Intensity slot have different number of samples"))

    # check RT slot
    if(nrt > 0) {
        if(nri != nrt)
            return(paste("Unequal number of RIs and RTs: ", nri,", ", nrt, sep = ""))
        nc.rt <- sapply(object@RT, ncol)
        if(!all(nc.rt == nc.rt[1]))
            return(paste("Some elements of the RT slot have different number of samples"))
    }
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
    if(any(dim(object@profRI)!=c(n,s)))
        paste("Unequal dimensions of 'profRI' slot.", dim(object@profRI)[1],", ", dim(object@profRI)[2], sep = "")
    if(any(dim(object@profInt)!=c(n,s)))
        paste("Unequal dimensions of 'profInt' slot.", dim(object@profInt)[1],", ", dim(object@profInt)[2], sep = "")
    if(nrow(object@info)!=n)
        paste("Unequal number of row of 'info' slot.", nrow(object@info), sep="")
    if(any(dim(object@profRI) > 0) & any(dim(object@profRI)!=c(n,s)))
        paste("Unequal dimensions of 'profRT' slot.", dim(object@profRT)[1],", ", dim(object@profRT)[2], sep = "")
     else  TRUE
})

    

setGeneric("as.list.tsProfile", function(x, ...) standardGeneric("as.list.tsProfile"))
setMethod("as.list.tsProfile", "tsMSdata", function(x, ...) as.list.tsMSdata(x))
