
# methods for tsMSdata
setGeneric("retIndex", function(obj) standardGeneric("retIndex"))
setMethod("retIndex", "tsMSdata", function(obj) obj@RI)
setGeneric("retIndex<-", function(obj, value) standardGeneric("retIndex<-"))
setReplaceMethod("retIndex", "tsMSdata", function(obj, value) { obj@RI <- value
 obj })

setGeneric("retTime", function(obj) standardGeneric("retTime"))
setMethod("retTime", "tsMSdata", function(obj) obj@RT)
setGeneric("retTime<-", function(obj, value) standardGeneric("retTime<-"))
setReplaceMethod("retTime", "tsMSdata", function(obj, value) { obj@RT <- value
 obj })

setGeneric("Intensity", function(obj) standardGeneric("Intensity"))
setMethod("Intensity", "tsMSdata", function(obj) obj@Intensity)
setGeneric("Intensity<-", function(obj, value) standardGeneric("Intensity<-"))
setReplaceMethod("Intensity", "tsMSdata", function(obj, value) { obj@Intensity <- value
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

# methods for tsProfile
setGeneric("profileInfo", function(obj) standardGeneric("profileInfo"))
setMethod("profileInfo", "tsProfile", function(obj) obj@info)
setGeneric("profileInfo<-", function(obj, value) standardGeneric("profileInfo<-"))
setReplaceMethod("profileInfo", "tsProfile", function(obj, value) { obj@info <- value
 obj })

setGeneric("profileInt", function(obj) standardGeneric("profileInt"))
setMethod("profileInt", "tsProfile", function(obj) obj@profInt)
setGeneric("profileInt<-", function(obj, value) standardGeneric("profileInt<-"))
setReplaceMethod("profileInt", "tsProfile", function(obj, value) { obj@profInt <- value
 obj })

setGeneric("profileRI", function(obj) standardGeneric("profileRI"))
setMethod("profileRI", "tsProfile", function(obj) obj@profRI)
setGeneric("profileRI<-", function(obj, value) standardGeneric("profileRI<-"))
setReplaceMethod("profileRI", "tsProfile", function(obj, value) { obj@profRI <- value
 obj })

setGeneric("profileRT", function(obj) standardGeneric("profileRT"))
setMethod("profileRT", "tsProfile", function(obj) obj@profRT)
setGeneric("profileRT<-", function(obj, value) standardGeneric("profileRT<-"))
setReplaceMethod("profileRT", "tsProfile", function(obj, value) { obj@profRT <- value
 obj })

setMethod("show", "tsProfile", function(object) {
	cat("An object of class 'tsMSdata':\n")
	cat(" Profile Information:\n")
	print(head(object@info), 5)
	if(nrow(object@info) > 5) cat("    ", nrow(object@info) - 5,"lines more...\n")
})

setGeneric("as.list.tsProfile", function(x, ...) standardGeneric("as.list.tsProfile"))
setMethod("as.list.tsProfile", "tsMSdata", function(x, ...) as.list.tsMSdata(x))
