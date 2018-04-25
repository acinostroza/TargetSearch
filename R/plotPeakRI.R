# select best peak according to method.
best_peak <- function(peaks, RI, method=c('RI', 'Intensity'))
{
    method <- match.arg(method)
    b <- split(1:nrow(peaks), peaks[,'fid'])
    best <- logical(nrow(peaks))
    if(method == 'Intensity') {
        a <- tapply(peaks[,'Int'], peaks[,'fid'], which.max)
    } else {
        a <- tapply(abs(RI - peaks[,'RI']), peaks[,'fid'], which.min)
    }
    # stopifnot(names(a) == names(b))
    for(i in 1:length(a))
        best[ b[[i]][ a[i] ] ] <- TRUE
    best
}

# function to approximate  RT/RI (y vector) from the RI/RT (x vector)
approxRT <- function(x, y) {
    m  <- coef(lm(y ~ x))
    ya <- pretty( m[1] + m[2]*x )
    xa <- (ya - m[1])/m[2]
    cbind(xa, ya)
}

`plotPeakRI` <-
function(smpInfo, refLib, id, dev=NULL, mz=NULL, RI=NULL,
         method=c('RI', 'Intensity'), useRI=TRUE, main=NULL,
         # graphical parameters:
         col=NULL, int_range=c(2,6), cex_range=c(.7,6), key_width=2)
{
    if(is_nullOrNA(RI))
        RI <- if(!is.na(medRI(refLib)[id])) medRI(refLib)[id] else libRI(refLib)[id]

    pk <- FindAllPeaks(smpInfo, refLib, id, dev, mz, RI, mz_type='quantMass')

    if(is.null(pk))
        return(invisible())
    if(nrow(pk) < 2)
        return(invisible())

    best <- best_peak(pk, RI, method)

    tmp <- rep(NA, length(smpInfo))
    names(tmp) <- sampleNames(smpInfo)
    tmp[pk[best, 'fid']] <- pk[best, 'RI']

    RTtime <- if(useRI) 'RT' else 'RI'
    RItime <- if(useRI) 'RI' else 'RT'

    mz   <- pk[1, 'mz']
    if(is.null(main)) {
        main <- sprintf("%s (%s) (%d)", libName(refLib)[id], names(libName(refLib)[id]), mz)
    }
    baseplot(pk[, 'fid'], pk[, RItime], log10(pk[, 'Int']), pk[, RTtime], best, main, RIexp=RI,
             nSamp=length(smpInfo), col=col, int_range=int_range, cex_range=cex_range, key_width=key_width)
    return(invisible(tmp))
}

# base plot
baseplot <- function(x, y, z, w, best, RIexp=NA, main="", nSamp=NULL,
                     col=NULL, int_range=c(2,6), cex_range=c(.7,6), key_width=2)
{
    if(is_nullOrNA(nSamp))
        nSamp <- max(x)
    if(is_nullOrNA(col))
        col <- c("#E0ECF4","#BFD3E6","#8C96C6","#88419D","#4D004B")
    rimed <- median(y[best])

    .scale <- function(x, a1, a2, b1, b2) (x - b1) / (b2 - b1) * (a2 - a1) + a1
    .colfun <- colorRamp(col, alpha=TRUE)
    .colrgb <- function(x) rgb(x[,1], x[,2], x[,3], x[,4], max=255)
    .limit  <- function(x, a, b) { x[x < a] <- a; x[x > b] <- b; x }

    zz   <- .limit(z, int_range[1], int_range[2])
    .cex <- .scale(zz, cex_range[1], cex_range[2], int_range[1], int_range[2])
    .bg  <- .colfun(.scale(zz, 0, 1, int_range[1], int_range[2]))
    .bg[!best, 4] <- 150
    .col <- .bg <- .colrgb(.bg)
    .pch <- sapply(best, function(x) if(x) 21 else 19)
    .col[best] <- 'black'

    layout(matrix(c(1,1,2,0), 2, 2), widths=c(1, lcm(key_width)))
    par(mar=c(5,4,4, 3)+.1)
    plot(x, y, cex=.cex, pch=.pch, bg=.bg, col=.col, ylab='RI', xlab='samples', main=main,
        panel.first=grid(), xlim=c(1,nSamp))
    abline(h=c(RIexp, rimed), lty=1:2)

    k <- setdiff(1:nSamp, x[best])
    if(length(k) > 0) {
        axis(1, at=k, labels=FALSE, tck=0.025)
        axis(3, at=k, labels=FALSE, tck=0.025)
    }

    tmp <- approxRT(y, w)
    axis(4, at=tmp[,1], label=tmp[,2])

    # gradient
    par(mar=c(0,0,4.1, 3.1))

    m    <- matrix(seq(0, 1, length=31), 1)
    mcol <- .colrgb(.colfun(m))
    image(m, col=mcol, axes=FALSE)
    box()

    pt <- pretty(int_range)
    axis(4, at = seq(0, 1, length.out=length(pt)), labels=format(10^pt, digits=2))
    invisible()
}

# vim: set ts=4 sw=4 et:
