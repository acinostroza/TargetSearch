# tsRim methods

setGeneric("rimMass", function(obj) standardGeneric("rimMass"))
setMethod("rimMass", "tsRim", function(obj) obj@mass)
setGeneric("rimMass<-", function(obj, value) standardGeneric("rimMass<-"))
setReplaceMethod("rimMass", "tsRim", function(obj, value) { obj@mass <- value
 validObject(obj)
 obj })

setGeneric("rimLimits", function(obj) standardGeneric("rimLimits"))
setMethod("rimLimits", "tsRim", function(obj) obj@limits)
setGeneric("rimLimits<-", function(obj, value) standardGeneric("rimLimits<-"))
setReplaceMethod("rimLimits", "tsRim", function(obj, value) { obj@limits <- value
 validObject(obj)
 obj })

setGeneric("rimStandard", function(obj) standardGeneric("rimStandard"))
setMethod("rimStandard", "tsRim", function(obj) obj@standard)
setGeneric("rimStandard<-", function(obj, value) standardGeneric("rimStandard<-"))
setReplaceMethod("rimStandard", "tsRim", function(obj, value) { obj@standard <- value
 validObject(obj)
 obj })

setValidity("tsRim", function(object) {
    if(ncol(object@limits) != 2)
        sprintf("Number of columns of limits (%d) it's not equal to 2.", ncol(object@limits))
    else if(nrow(object@limits) != length(object@standard))
        sprintf("Unequal number of standards (%d) and limits (%d)",
                length(object@standard), nrow(object@limits))
    else if(length(object@mass) != 1 && length(object@mass) != length(object@standard))
        sprintf("Unequal number of standards (%d) and mass markers (%d)",
                length(object@standard), length(object@mass))
    else TRUE
})

setMethod("[", "tsRim",
    function(x, i, j, ..., drop=TRUE)
    {
        if(missing(i) & missing(j)) {
            validObject(x)
            return(x)
        }
        if(!missing(j))
            stop("incorrect number of dimensions")

        if(is.character(i)) {
            if(is.null(rownames(x@limits)))
                stop("subscript out of bounds")
            i <- match(i, rownames(x@limits))
        }
        assert_that(noNA(i), msg="subscript out of bounds")

        x@limits <- x@limits[i,,drop=FALSE]
        x@standard <- x@standard[i]
        x@mass <- if(length(x@mass) == 1) x@mass else x@mass[i]
        validObject(x)
        return(x)
    }
)

setMethod("c", "tsRim",
    function(x, ..., recursive=FALSE)
    {
        assert_that(is.flag(recursive), !recursive,
                    msg="\"c\" method for `tsRim` objects does not support the 'recursive' option")
        rep2 <- function(x, n)
            if(length(x) == 1) rep(x, n) else x

        z <- list(...)
        if(length(z) == 0)
            return(x)

        z <- c(list(x), z)
        lim <- lapply(z, slot, 'limits')
        n   <- vapply(lim, nrow, 0)
        lim <- do.call('rbind', lim)
        std <- unlist(lapply(z, slot, 'standard'))
        mass <- lapply(z, slot, 'mass')
        m <- unique(unlist(mass))
        mass <- if(length(m) == 1) m else
                  unlist(mapply(rep2, mass, n, SIMPLIFY=FALSE))
        obj <- new('tsRim', limits=lim, standard=std, mass=mass)
        validObject(obj)
        return(obj)
    }
)

setMethod("length", "tsRim", function(x) length(x@standard))

setMethod("initialize", "tsRim", function(.Object, limits, standard, mass, ...)
{
    if(!is.numeric(limits))
        stop("The RI marker limit definition must be a numeric matrix")

    if(is.null(rownames(limits)) && !is.null(dim(limits)))
        rownames(limits) <- paste("RI.Marker", seq(nrow(limits)))

    .Object@limits <- limits
    .Object@standard <- standard
    .Object@mass <- mass
    validObject(.Object)
    .Object
})

# vim: set ts=4 sw=4 et:
