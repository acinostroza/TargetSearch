`quantMatrix` <-
    function(Lib, metabProfile, value=c("quantmass", "maxint", "maxobs"),
             selmass=FALSE)
{
    value <- match.arg(value)
    assert_that(is.flag(selmass))
    Int   <- Intensity(metabProfile)
    id <- rownames(profileInfo(metabProfile))
    stopifnot(id %in% names(selMass(Lib)))

    # match library and profile unique identifiers
    Lib <- Lib[ id ]

    # sanity checks. rownames of profiles must be equal to top masses
    # of library
    rnames <- lapply(Int, rownames)
    equal <- function(a, b) length(a) == length(b) && all(as.character(a) == b)
    tmp <- mapply(equal, topMass(Lib), rnames)
    if(!all(tmp))
        stop("Library and Profile objects do not match")

    # quant mass must be first row
    Q <- vapply(Int, function(x) as.numeric(rownames(x)[1]), 0)
    if(!all(Q == quantMass(Lib)))
        stop("Library and Profile objects do not match")

    # if selmass is TRUE, restrict rows to selective masses only
    if(selmass)
        Int <- mapply(function(x, m) x[ as.character(m), ],
                      Int, selMass(Lib), SIMPLIFY=FALSE)

    sM <- strsplit(profileInfo(metabProfile)$Masses, ";")
    sM <- lapply(sM, as.numeric)

    if(value == 'quantmass') {
        best <- rep.int(1, length(Int))
    }

    else if(value == "maxint") {
        best <- vapply(Int, function(x) {
            z <- apply(x, 1, median, na.rm=TRUE)
            which.max(z)
        }, 0)
    }

    else if(value == "maxobs") {
        best <- vapply(Int, function(x) {
            z <- rowSums(!is.na(x))
            which.max(z)
        }, 0)
    }

    else {
        stop("Invalid parameter `value`: ", value)
    }

    M <- t(mapply(function(x, k) x[k,], Int, best))
    Q <- mapply(function(x, k) rownames(x)[k], Int, best)
    Q <- as.numeric(Q)
    f <- function(a, b) a %in% b
    attr(M, "quantMass") <- Q
    attr(M, "isSelMass") <- mapply(f, Q, selMass(Lib))
    attr(M, "isCorMass") <- mapply(f, Q, sM)
    attr(M, "libNames") <- profileInfo(metabProfile)$Name
    M
}

# vim: set ts=4 sw=4 et:
