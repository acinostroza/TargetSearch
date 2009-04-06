
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

setMethod("show", "tsMSdata", function(object) {
	cat("An object of class 'tsMSdata':\n")
	cat(" Number of samples: ", ncol(object@RI), "\n")
	cat(" Number of masses:  ", nrow(object@RI), "\n")
})

# methods for tsProfile
setGeneric("profileInfo", function(obj) standardGeneric("profileInfo"))
setMethod("profileInfo", "tsProfile", function(obj) obj@info)
setGeneric("profileInfo<-", function(obj, value) standardGeneric("profileInfo<-"))
setReplaceMethod("profileInfo", "tsProfile", function(obj, value) { obj@info <- value
 obj })

setMethod("show", "tsProfile", function(object) {
	cat("An object of class 'tsMSdata':\n")
	cat(" Profile Information:\n")
	print(head(object@info), 5)
	if(nrow(object@info) > 5) cat("    ", nrow(object@info) - 5,"lines more...\n")
})

