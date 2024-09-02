
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
 obj <- .setLibrary(obj)
 obj })

setMethod("length", "tsLib", function(x) length(x@medRI))

setGeneric("selMass", function(obj) standardGeneric("selMass"))
setMethod("selMass", "tsLib", function(obj) obj@selMass)
setGeneric("selMass<-", function(obj, value) standardGeneric("selMass<-"))
setReplaceMethod("selMass", "tsLib", function(obj, value) { obj@selMass <- value
 obj <- .setLibrary(obj)
 obj })


setGeneric("topMass", function(obj) standardGeneric("topMass"))
setMethod("topMass", "tsLib", function(obj) obj@topMass)
setGeneric("topMass<-", function(obj, value) standardGeneric("topMass<-"))
setReplaceMethod("topMass", "tsLib", function(obj, value) { obj@topMass <- value
 obj <- .setLibrary(obj)
 obj })

setGeneric("quantMass", function(obj) standardGeneric("quantMass"))
setMethod("quantMass", "tsLib", function(obj) obj@quantMass)
setGeneric("quantMass<-", function(obj, value) standardGeneric("quantMass<-"))
setReplaceMethod("quantMass", "tsLib", function(obj, value) { obj@quantMass <- as.numeric(value)
 obj <- .setLibrary(obj)
 obj })


setGeneric("spectra", function(obj) standardGeneric("spectra"))
setMethod("spectra", "tsLib", function(obj) obj@spectra)
setGeneric("spectra<-", function(obj, value) standardGeneric("spectra<-"))
setReplaceMethod("spectra", "tsLib", function(obj, value) {
    if(is.null(value) || length(value) == 0)
        obj@spectra <- vector("list", length(obj))
    else if(is.list(value) && length(obj) == length(value))
        obj@spectra <- value
    else
        stop("Invalid value. `spectra` must be a ", length(obj), "-list of matrices.")
    .setLibrary(obj)
})

setGeneric("libName", function(obj) standardGeneric("libName"))
setMethod("libName", "tsLib", function(obj) obj@Name)
setGeneric("libName<-", function(obj, value) standardGeneric("libName<-"))
setReplaceMethod("libName", "tsLib", function(obj, value) { obj@Name <- value
 obj <- .setLibrary(obj)
 obj })

setGeneric("libRI", function(obj) standardGeneric("libRI"))
setMethod("libRI", "tsLib", function(obj) obj@RI)
setGeneric("libRI<-", function(obj, value) standardGeneric("libRI<-"))
setReplaceMethod("libRI", "tsLib", function(obj, value) { obj@RI <- value
 obj <- .setLibrary(obj)
 obj })

setGeneric("libData", function(obj) standardGeneric("libData"))
setMethod("libData", "tsLib", function(obj) obj@libData)
setGeneric("libData<-", function(obj, value) standardGeneric("libData<-"))
setReplaceMethod("libData", "tsLib", function(obj, value) {
	obj@libData <- value
	obj <- .setLibrary(obj)
	obj
})


setGeneric("RIdev", function(obj) standardGeneric("RIdev"))
setMethod("RIdev", "tsLib", function(obj) obj@RIdev)
setGeneric("RIdev<-", function(obj, value) standardGeneric("RIdev<-"))
setReplaceMethod("RIdev", "tsLib", function(obj, value) {
 obj@RIdev <- value
 obj <- .setLibrary(obj)
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
    x@RIdev <- x@RIdev[i,,drop=FALSE]
    x@selMass <- x@selMass[i]
    x@topMass <- x@topMass[i]
    x@libData <- x@libData[i, j, drop=FALSE]
    x@spectra <- x@spectra[i]
    x@quantMass <- x@quantMass[i]
    .setLibrary(x)
})

setValidity("tsLib", function(object) {
    check_spectra <- function(spectra) {
        check <- function(x)
            (is.null(x)) || (length(x) == 0) || (is.numeric(x) && is.matrix(x) && (ncol(x) == 2))
        all(sapply(spectra, check))
    }
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
	else if(length(object@spectra) != n)
		paste("Unequal number of Names and spectra: ", n,", ", length(object@spectra), sep = "")
	else if(nrow(object@libData) != n)
		paste("Unequal number of Names and libData: ", n,", ", nrow(object@libData), sep = "")
    else if(!check_spectra(object@spectra))
        paste("Detected invalid elements in `spectra`. Expecting a list of two-column matrices.")
	else TRUE
})

setMethod("$", "tsLib", function(x, name) {
    eval(substitute(libData(x)$NAME_ARG, list(NAME_ARG=name)))
})

setMethod("c", "tsLib",
    function(x, ..., recursive=FALSE)
    {
        assert_that(is.flag(recursive), !recursive,
                    msg="\"c\" method for `tsLib` objects does not support the 'recursive' option")

        z <- list(...)
        if(length(z) == 0)
            return(x)

        z <- c(list(x), z)
        stopifnot(all(vapply(z, validObject, FALSE)))
        Name <- unlist(lapply(z, slot, 'Name'))
        RI <- unlist(lapply(z, slot, 'RI'))
        medRI <- unlist(lapply(z, slot, 'medRI'))
        RIdev <- do.call('rbind', lapply(z, slot, 'RIdev'))
        selMass <- unlist(lapply(z, slot, 'selMass'), recursive=FALSE)
        topMass <- unlist(lapply(z, slot, 'topMass'), recursive=FALSE)
        quantMass <- unlist(lapply(z, slot, 'quantMass'))
        spectra <- unlist(lapply(z, slot, 'spectra'), recursive=FALSE)

        # try first to merge libData, otherwiwse fails with a more meaningful message
        libData <- try(do.call(.rbind, lapply(z, slot, 'libData')), silent=TRUE)
        if(inherits(libData, 'try-error'))
            stop('Unable to combine `libData` slot due to incompatible data types. ',
                 'Please call the `libData` method on the objects to combine and check ',
                 'they are compatible.')

        obj <- new('tsLib', Name=Name, RI=RI, selMass, medRI=medRI, RIdev=RIdev,
                   topMass=topMass, quantMass=quantMass, spectra=spectra,
                   libData=libData)
        validObject(obj)
        obj
    }
)

setMethod("initialize",
          "tsLib",
          function(.Object, Name, RI, selMass, medRI=RI, RIdev=NULL,
                            topMass=NULL, quantMass=NULL, spectra=vector("list", length(Name)),
                            libData=NULL)
          {
            # require at least 1 Name, 1 RI and 1 selMass
            if (length(Name) != length(RI) | length(Name) != length(selMass))
                stop("'Name', 'RI', and 'selMass' must have the same length")

            if(is.null(RIdev)) {
                RIdev <- matrix(rep(RI, 3), ncol=3)
                RIdev <- sweep(matrix(rep(RI, 3), ncol=3), 2, c(5,2,1)/1000, FUN="*")
            } else if(!is.numeric(RIdev)) {
                stop("'RIdev' must be a numeric 3-column matrix")
            } else if(!is.matrix(RIdev)) {
                if(length(RIdev) == 3) {
                    RIdev <- matrix(RIdev, nrow=length(Name), ncol=3, byrow=TRUE)
                } else if(length(RIdev) == 3 * length(Name)) {
                    RIdev <- matrix(RIdev, nrow=length(Name), ncol=3)
                }
            }

            if(is.numeric(selMass)) {
                selMass <- as.list(selMass)
            }

            if(!is.numeric(sapply(selMass, getElement, 1))) {
                stop("'selMass' must be a list of numeric vectors")
            }

            if(is.null(topMass))
                topMass <- selMass
            if(is.null(quantMass))
                quantMass <- sapply(selMass, getElement, 1)
            if(is.null(libData))
                libData <- data.frame(Name = Name, RI = RI)

            if(is.null(spectra))
                spectra <- vector("list", length(Name))

            .Object@Name     <- Name
            .Object@RI       <- RI
            .Object@medRI    <- medRI
            .Object@RIdev    <- RIdev
            .Object@selMass  <- selMass
            .Object@topMass  <- topMass
            .Object@quantMass <- quantMass
            .Object@spectra  <- spectra
            .Object@libData  <- libData
            .Object <- .setLibrary(.Object)
            .Object
          })

.setLibrary <- function(lib)
{
    # remove NAs and take unique
    uniqrm <- function(x) {
        y <- x[ !is.na(x) ]
        unique(y)
    }

    ma <- function(...) mapply(..., SIMPLIFY=FALSE)
    uf <- function(...) uniqrm(c(...))

    id <- lib@libData$libID

    if(is.null(id)) {
        lib@libData <- .addLibID(lib@libData)
        id <- lib@libData$libID
    }

    # remove names in data.frame libData
    dat <- lapply(as.list( lib@libData ), function(x) { unname(x) })
    dat <- data.frame(dat, stringsAsFactors=FALSE)

    lib@selMass <- ma(uf, lib@quantMass, lib@selMass)
    lib@topMass <- ma(uf, lib@quantMass, lib@selMass, lib@topMass)
    lib@quantMass <- sapply(lib@topMass, getElement, 1)

    rownames(lib@RIdev) <- rownames(dat) <- id
    colnames(lib@RIdev) <- sprintf("Win_%d", 1:3)
    names(lib@Name) <- names(lib@RI) <- names(lib@medRI) <- names(lib@selMass) <- id
    names(lib@topMass) <- names(lib@quantMass) <- names(lib@spectra) <- id

    # remove extra columns from dat
    extra <- c("Win_1", "Win_2", "Win_3", "SEL_MASS", "TOP_MASS", "SPECTRUM", "QUANT_MASS")
    k <- setdiff(colnames(dat), extra)
    dat <- dat[, k,drop=FALSE]

    lib@libData <- dat
    stopifnot( validObject(lib) )
    lib
}

.addLibID <- function(lib) {
    k <- which(tolower(colnames(lib)) == "libid")
    if(length(k) == 0) {
        lib <- data.frame(libID=paste("GC", 1:nrow(lib), sep="."), lib, stringsAsFactors=FALSE)
    } else if(length(k) == 1) {
        if(colnames(lib)[k] != 'libID') {
            warning(sprintf("Changing '%s' to 'libID'", colnames(lib)[k]))
            colnames(lib)[k] <- 'libID'
        }
    }
    else {
        stop("Multiple colnames match 'libID'. Expecting exactly one match. ",
             "Please rename/remove the extra columns")
    }
    id <- make.names(lib$libID, TRUE)
    if(any(id != lib$libID))
        warning("Some identifiers where renamed in order to make them unique")
    lib$libID <- id
    return(lib)
}

# vim: set ts=4 sw=4 et:
