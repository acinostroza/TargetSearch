findRetentionTime <- function(rTime, Intensity, limits)
{
    low <- findInterval(limits[,1], rTime) + 1
    upp <- findInterval(limits[,2], rTime)
    get <- function(x, k) if(!is.null(dim(x))) x[,k] else x

    ret <- mapply(function(a, b, k) {
                    if(a <= b) (which.max(get(Intensity, k)[a:b]) + a - 1) else NA
                  }, low, upp, seq(nrow(limits)))
    int <- mapply(function(i, j) get(Intensity, j)[i], ret, seq(nrow(limits)))
    structure(rTime[ ret ], intensity=Intensity[ ret ])
}
