# tsRim methods

setGeneric("rimMass", function(obj) standardGeneric("rimMass"))
setMethod("rimMass", "tsRim", function(obj) obj@mass)
setGeneric("rimMass<-", function(obj, value) standardGeneric("rimMass<-"))
setReplaceMethod("rimMass", "tsRim", function(obj, value) { obj@mass <- value
 obj })
 
setGeneric("rimLimits", function(obj) standardGeneric("rimLimits"))
setMethod("rimLimits", "tsRim", function(obj) obj@limits)
setGeneric("rimLimits<-", function(obj, value) standardGeneric("rimLimits<-"))
setReplaceMethod("rimLimits", "tsRim", function(obj, value) { obj@limits <- value
 obj })

setGeneric("rimStandard", function(obj) standardGeneric("rimStandard"))
setMethod("rimStandard", "tsRim", function(obj) obj@standard)
setGeneric("rimStandard<-", function(obj, value) standardGeneric("rimStandard<-"))
setReplaceMethod("rimStandard", "tsRim", function(obj, value) { obj@standard <- value
 obj })

setValidity("tsRim", function(object) {
	if(ncol(object@limits) != 2)
		message("Number of columns of limits it is not equal to 2")
	else if(nrow(object@limits) != length(object@standard)) 
		paste("Unequal number of standards and limits: ", length(object@standard),
			", ", nrow(object@limits), sep = "")
	else TRUE
})


