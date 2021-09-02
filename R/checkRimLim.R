#' Visually check retention index marker limits
#'
#' A function to visually check if the retention time limits of the retention
#' index markers (aka FAMEs) are correct.
#'
#' @param samples A  \code{\linkS4class{tsSample}} object created by
#'    \code{\link{ImportSamples}}.
#' @param rim A \code{\linkS4class{tsRim}} object describing the retention index
#'    markers. See \code{\link{ImportFameSettings}}.
#' @param layout A vector of the form \code{c(nr, nc)} to arrange the panel by
#'    \code{nr} rows and \code{nc} columns. If missing then the layout is created
#'    automatically.
#' @param single Logical. If TRUE, a single sample will be selected randomly
#'    for plotting. This was the old default behavior. If FALSE, all samples
#'    will be used for plotting (see note below).
#' @param show Logical. If \code{FALSE} the plot is not shown, but the data points
#'    can be used for further inspection or for custom plots.
#' @param extend a numeric coefficient to extend the time window search of the
#'    respective time marker. Defaults to \code{0.5}.
#' @param rect.col the color for the background rectangle which indicates
#'    the current retention time limits.
#' @param mar the subplots margins, passed to \code{\link[graphics]{par}()}.
#' @param oma the outer plot margins, passed to \code{\link[graphics]{par}()}.
#' @param cex.main The magnification to be used for main titles, passed to
#'    \code{\link[graphics]{par}()}.
#' @param type A character vector indicating the type of plots. Default \code{"l"} for
#'    lines. Passed to \code{\link[base]{plot}()}.
#' @param ... extra plotting arguments passed to \code{\link[base]{plot}()}
#'    such as \code{col}, \code{lty}, \code{pch}, \code{lwd}.
#'
#' @return
#'    A list of \code{n} times \code{2} matrices or invisible. Each element correspond to a
#'    marker. Columns are retention time and intensities of the respective marker's m/z.
#'    The rows can be as many data points are within the search window.
#'
#' @seealso
#'    \code{\linkS4class{tsSample}}, \code{\linkS4class{tsRim}}, \code{\link{ImportFameSettings}}
`checkRimLim` <-
function(samples, rim, layout, show=TRUE, single=TRUE, extend=0.5,
    rect.col="#e7e7e7", mar=c(2,2,2,2), oma=c(3,3,2,0.5), cex.main=1, type='l', ...)
{
    panel <- function(z, r) {
        ylim <- range(z)
        ylim <- ylim + c(-1, 1)*0.04*diff(ylim)
        rect(r[1], ylim[1], r[2], ylim[2], col=rect.col, border=FALSE)
    }
    # single element
    se <- function(x) if(length(x) == 1) x[[1]] else x
    assert_that(is.tsSample(samples))
    assert_that(is.tsRim(rim))
    assert_that(is.flag(show))
    assert_that(is.flag(single))

    if(single)
        samples <- samples[ sample(length(samples), 1) ]
    nr <- nrow(rim@limits)
    rim@mass <- if(length(rim@mass) == 1) rep(rim@mass, nr) else rim@mass

    if(missing(layout)) {
        n  <- ceiling(sqrt(nr))
        layout <- if(n * (n - 1) >= nr) c(n-1, n) else c(n, n)
    } else {
        assert_that(is.numeric(layout), length(layout) == 2)
        if(prod(layout) < nr)
            stop(sprintf("Invalid layout. The supplied value cannot hold %d panels", nr))
    }

    dat <- lapply(CDFfiles(samples), .chkrl_get_data, rim, extend)
    names(dat) <- sampleNames(samples)

    if(!show)
	    return(invisible(se(dat)))

    pdata <- .chkrl_prep_data(dat, nr)
    op <- par(mfrow = layout, mar=mar, oma=oma, cex.main=cex.main)
    on.exit(par(op))
    for(i in 1:nr) {
        Time <- pdata[[i]][, 1] ; Intensity <- pdata[[i]][, -1]
        matplot(Time, Intensity, type=type, panel.first=panel(Intensity, rim@limits[i, ]),
            main=rownames(rim@limits)[i], ...)
        legend('topright', legend=sprintf("m/z: %d", rim@mass[i]), box.lty=0)
    }
    title(main = .chkrl_main(names(dat)), outer=TRUE)

    invisible(se(dat))
}

# checkRimLim Helper functions
# all functions are prefixed with .chkrl as they are used only here

.chkrl_get_data <- function(cdf, ...)
{
    if(.is_ts_ncdf4(cdf)) .chkrl_get_data4(cdf, ...) else .chkrl_get_data3(cdf, ...)
}

# get data from netcdf4 files
.chkrl_get_data4 <- function(cdf, rim, extend=0.5)
{
    nc <- nc_open(cdf)
    on.exit(nc_close(nc))
    Time <- ncvar_get(nc, 'retention_time')
    mzRange <- ncvar_get(nc, 'mass_range')
    ret <- mapply(function(low, up, mz) {
        y <- c(low, up) + (up - low) * c(-1, 1) * extend
        z <- which(Time > y[1] & Time < y[2])
        w <- ncvar_get(nc, 'intensity', start=c(z[1], mz - mzRange[1] + 1), count=c(length(z), 1))
        cbind(Time=Time[z], Intensity=w)
        }, rim@limits[,1], rim@limits[,2], rim@mass, SIMPLIFY=FALSE)
}

# get data from netcdf3 files
.chkrl_get_data3 <- function(cdf, rim, extend=0.5)
{
    dat <- peakCDFextraction(cdf)
    mzRange <- dat$massRange
    Time <- dat$Time
    ret <- mapply(function(low, up, mz, extend) {
        y <- c(low, up) + (up - low) * c(-1, 1) * extend
        z <- which(Time > y[1] & Time < y[2])
        w <- dat$Peaks[ z, mz - mzRange[1] + 1 ]
        cbind(Time=Time[z], Intensity=w)
        }, rim@limits[,1], rim@limits[,2], rim@mass, extend, SIMPLIFY=FALSE)
}

# prepare data for multiple plotting
.chkrl_prep_data <- function(dat, nr) {
    res <- lapply(seq(nr), function(j) {
        z <- lapply(dat, getElement, j)
        t <- range(sapply(z, function(x) range(x[,1])))
        n <- max(sapply(z, nrow))
        t <- seq(t[1], t[2], length=n)
        y <- sapply(z, function(x) approx(x[,1], x[,2], t)$y)
        cbind(time=t, y)
    })
    res
}

# make title for plot
.chkrl_main <- function(x, len=6) {
    s <- if(length(x) == 1) "" else "s"
    if(length(x) > len)
        paste0("Samples: ", paste(x[1:len], collapse=", "), ", ...")
    else
        paste0("Sample", s, ": ", paste(x, collapse=", "))
}
