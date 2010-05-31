
# tsLib	
setGeneric("refLib", function(obj, ...) standardGeneric("refLib"))
setMethod("refLib", "tsLib", function(obj, ri, w = 1, sel = TRUE) {
	if(missing(ri))
		ri <- obj@medRI
	
	if(sel)	{
		rp <- sapply(obj@selMass,length)
		m  <- unlist(obj@selMass)
	} else {
		rp <- sapply(obj@topMass,length)
		m  <- unlist(obj@topMass)
	}
	win <- rep(obj@RIdev[,w], rp)
	ri2 <- rep(ri, rp)
	out <- cbind(ri2 - win, m, ri2 + win)
	colnames(out) <- c("minRI", "mz", "maxRI")
	rownames(out) <- libId(obj, sel)
	out
})

setGeneric("libId", function(obj, ...) standardGeneric("libId"))
setMethod("libId", "tsLib", function(obj, sel = TRUE) {
	if(sel)
		rep(1:length(obj@RI), sapply(obj@selMass, length))
	else
		rep(1:length(obj@RI), sapply(obj@topMass, length))
})

setGeneric("medRI", function(obj) standardGeneric("medRI"))
setMethod("medRI", "tsLib", function(obj) obj@medRI)
setGeneric("medRI<-", function(obj, value) standardGeneric("medRI<-"))
setReplaceMethod("medRI", "tsLib", function(obj, value) { obj@medRI <- value
 validObject(obj)
 obj })

setMethod("length", "tsLib", function(x) length(x@medRI))

setGeneric("selMass", function(obj) standardGeneric("selMass"))
setMethod("selMass", "tsLib", function(obj) obj@selMass)
setGeneric("selMass<-", function(obj, value) standardGeneric("selMass<-"))
setReplaceMethod("selMass", "tsLib", function(obj, value) { obj@selMass <- value
 validObject(obj)
 obj })


setGeneric("topMass", function(obj) standardGeneric("topMass"))
setMethod("topMass", "tsLib", function(obj) obj@topMass)
setGeneric("topMass<-", function(obj, value) standardGeneric("topMass<-"))
setReplaceMethod("topMass", "tsLib", function(obj, value) { obj@topMass <- value
 validObject(obj)
 obj })

setGeneric("quantMass", function(obj) standardGeneric("quantMass"))
setMethod("quantMass", "tsLib", function(obj) obj@quantMass)
setGeneric("quantMass<-", function(obj, value) standardGeneric("quantMass<-"))
setReplaceMethod("quantMass", "tsLib", function(obj, value) { obj@quantMass <- as.numeric(value)
 validObject(obj)
 obj })


setGeneric("spectra", function(obj) standardGeneric("spectra"))
setMethod("spectra", "tsLib", function(obj) obj@spectra)
setGeneric("spectra<-", function(obj, value) standardGeneric("spectra<-"))
setReplaceMethod("spectra", "tsLib", function(obj, value) { obj@spectra <- value
 validObject(obj)
 obj })

setGeneric("libName", function(obj) standardGeneric("libName"))
setMethod("libName", "tsLib", function(obj) obj@Name)
setGeneric("libName<-", function(obj, value) standardGeneric("libName<-"))
setReplaceMethod("libName", "tsLib", function(obj, value) { obj@Name <- value
 validObject(obj)
 obj })

setGeneric("libRI", function(obj) standardGeneric("libRI"))
setMethod("libRI", "tsLib", function(obj) obj@RI)
setGeneric("libRI<-", function(obj, value) standardGeneric("libRI<-"))
setReplaceMethod("libRI", "tsLib", function(obj, value) { obj@RI <- value
 validObject(obj)
 obj })

setGeneric("libData", function(obj) standardGeneric("libData"))
setMethod("libData", "tsLib", function(obj) obj@libData)
setGeneric("libData<-", function(obj, value) standardGeneric("libData<-"))
setReplaceMethod("libData", "tsLib", function(obj, value) {
	obj@libData <- value
	validObject(obj)
	obj
})


setGeneric("RIdev", function(obj) standardGeneric("RIdev"))
setMethod("RIdev", "tsLib", function(obj) obj@RIdev)
setGeneric("RIdev<-", function(obj, value) standardGeneric("RIdev<-"))
setReplaceMethod("RIdev", "tsLib", function(obj, value) {
 obj@RIdev <- value
 validObject(obj)
 obj
})

setMethod("show", "tsLib", function(object) {
	cat("An object of class 'tsLib':\n")
	cat(" Number of objects:  ", length(object), "\n")
	cat("\nImported Library Info:\n")
	print(head(libData(object), 5))
	if(length(object) > 5) cat("    ", length(object) - 5,"lines more...\n")
})

setMethod("[", "tsLib", function(x, i, j, ..., drop) {
    x@Name <- x@Name[i];
    x@RI <- x@RI[i];
    x@medRI <- x@medRI[i];
    x@RIdev <- matrix(x@RIdev[i,], ncol = 3)
    x@selMass <- x@selMass[i]
    x@topMass <- x@topMass[i]
    x@libData <- x@libData[i,]
    x@spectra <- x@spectra[i]
    x@quantMass <- x@quantMass[i]
    x
})

setValidity("tsLib", function(object) {
	n <- length(object@Name)
	if(length(object@RI) != n)
		paste("Unequal number of Names and RI: ", n,", ", length(object@RI), sep = "")
	else if(length(object@medRI) != n)
		paste("Unequal number of Names and medRI: ", n,", ", length(object@medRI), sep = "")
	else if(nrow(object@RIdev) != n)
		paste("Unequal number of Names and RIdev: ", n,", ", nrow(object@RIdev), sep = "")
	else if(ncol(object@RIdev) != 3)
		paste("Number of columns of RIdev is not 3: ", ncol(object@RIdev), sep = "")
	else if(length(object@selMass) != n)
		paste("Unequal number of Names and selMass: ", n,", ", length(object@selMass), sep = "")
	else if(length(object@topMass) != n)
		paste("Unequal number of Names and topMass: ", n,", ", length(object@topMass), sep = "")
	else if(length(object@quantMass) != n)
		paste("Unequal number of Names and quantMass: ", n,", ", length(object@quantMass), sep = "")
	else if(length(object@spectra) != n & length(object@spectra) != 0)
		paste("Unequal number of Names and spectra: ", n,", ", length(object@spectra), sep = "")
	else if(nrow(object@libData) != n)
		paste("Unequal number of Names and libData: ", n,", ", nrow(object@libData), sep = "")
	else TRUE
})

setMethod("$", "tsLib", function(x, name) {
    eval(substitute(libData(x)$NAME_ARG, list(NAME_ARG=name)))
})

setMethod("initialize",
          "tsLib",
          function(.Object, Name    = character(0), RI      = numeric(0),
                            medRI   = numeric(0),   RIdev   = matrix(0,0,3),
                            selMass = list(),       topMass = list(),
                            quantMass = numeric(0),
                            spectra = list(),       libData = data.frame()) {
            if (length(Name) > 0) {
                if(length(medRI) == 0)
                    medRI   <- RI
                if(length(topMass) == 0)
                    topMass <- selMass
                if(length(quantMass) == 0)
                    quantMass <- numeric(length(Name))
                if(all(dim(libData) == 0))
                    libData <- data.frame(Name = Name, RI = RI)
                if(is.matrix(RIdev) == FALSE & length(RIdev) == 3)
                    RIdev   <- matrix(rep(RIdev, length(Name)), length(Name), 3, byrow = T)
            }
            ids <- as.character(1:length(Name))
            .Object@Name     <- Name
            .Object@RI       <- RI
            .Object@medRI    <- medRI
            .Object@RIdev    <- RIdev
            .Object@selMass  <- selMass
            .Object@topMass  <- topMass
            .Object@quantMass <- quantMass
            .Object@spectra  <- spectra
            .Object@libData  <- libData
            names(.Object@Name)     <- ids
            names(.Object@RI)       <- ids
            names(.Object@medRI)    <- ids
            rownames(.Object@RIdev) <- ids
            names(.Object@selMass)  <- ids
            names(.Object@topMass)  <- ids
            if(length(.Object@spectra) > 0) names(.Object@spectra)  <- ids
            if(length(.Object@quantMass) > 0) names(.Object@quantMass)  <- ids
            .Object
          })
