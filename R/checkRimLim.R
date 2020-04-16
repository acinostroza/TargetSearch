#' Visually check retention index marker limits
#'
#' A function to visually check if the retention time limites of the retention
#' index markers (aka FAMEs) is correct.
#'
#' The function takes a random CDF file from your \code{\linkS4class{tsSample}}
#' object and creates a panel plot of the m/z traces around the area in which a
#' marker is expected to be. Repeated calls to this function can be used to check
#' other samples. It is also possible to check a specific sample by indexing
#' the \code{\linkS4class{tsSample}} object.
#'
#' @param samples A  \code{\linkS4class{tsSample}} object created by
#'    \code{\link{ImportSamples}}.
#' @param rim A \code{\linkS4class{tsRim}} object describing the retention index
#'    markers. See \code{\link{ImportFameSettings}}.
#' @param layout A vector of the form \code{c(nr, nc)} to arrange the panel by
#'    \code{nr} rows and \code{nc} columns. If missing then the layout is created
#'    automatically.
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
#'    lines. Passed to \code{\link[graphics]{plot}()}.
#' @param ... extra plotting arguments passed to \code{\link[graphics]{plot}()}
#'    such as \code{col}, \code{lty}, \code{pch}, \code{lwd}.
#'
#' @return
#'    A list of \code{n} times \code{2} matrices or invisible. Each element correspond to a
#'    marker. Columns are retention time and intensities of the respective marker's m/z.
#'    The rows can be as many data points are within the search window.
#'
#' @seealso
#'    \code{\linkS4class{tsSample}}, \code{\linkS4class{tsRim}}, \code{\link{ImportFameSettings}}
#'
#' @examples
#' require(TargetSearchData)
#'
#' # get the cdf path TargetSearchData
#' cdfpath <- file.path(find.package("TargetSearchData"), "gc-ms-data")
#'
#' # import samples (see ImportSamples() for details)
#' samples <- ImportSamples(file.path(cdfpath, "samples.txt"), CDFpath = cdfpath)
#'
#' # Import RI markers (see ImportFameSettings())
#' rim <- ImportFameSettings(file.path(cdfpath, "rimLimits.txt"))
#'
#' # choose a sample at random and plot the m/z traces around the retention
#' # time window
#' ret <- checkRimLim(samples, rim)
#'
#' # to choose a specific samples and marker, use subsetting
#' ret <- checkRimLim(samples[3], rim[2])
#'
`checkRimLim` <-
function(samples, rim, layout, show=TRUE, extend=0.5, rect.col="#e7e7e7",
    mar=c(2,2,2,2), oma=c(3,3,2,0.5), cex.main=1, type='l', ...)
{
    panel <- function(z, r) {
        ylim <- range(z[, 2])
        ylim <- ylim + c(-1, 1)*0.04*diff(ylim)
        rect(r[1], ylim[1], r[2], ylim[2], col=rect.col, border=FALSE)
    }
    k <- sample(length(samples), 1)
    nr <- nrow(rim@limits)
    rim@mass <- if(length(rim@mass) == 1) rep(rim@mass, nr) else rim@mass

    if(missing(layout)) {
        n  <- ceiling(sqrt(nr))
        layout <- if(n * (n - 1) >= nr) c(n-1, n) else c(n, n)
    } else {
        layout <- layout[1:2]
        if(prod(layout) < nr)
            stop("invalid layout")
    }

    cdf <- CDFfiles(samples)[k]
    get_data <- if(.is_ts_ncdf4(cdf)) .get_data4 else .get_data3
    dat <- get_data(cdf, rim, extend)
    names(dat) <- rownames(rim@limits)

    if(!show)
	    return(invisible(dat))

    mass <- rim@mass
    if(length(mass) == 1)
        mass <- rep(mass, nrow(rim@limits))

    op <- par(mfrow = layout, mar=mar, oma=oma, cex.main=cex.main)
    for(i in 1:nr) {
        plot(dat[[i]], type=type, panel.first=panel(dat[[i]],rim@limits[i, ]),
            main=rownames(rim@limits)[i], ...)
        legend('topright', legend=sprintf("m/z: %d", mass[i]), box.lty=0)
    }
    title(main = basename(cdf), outer=TRUE)
    par(op)
    invisible(dat)
}

# get data from netcdf4 files
.get_data4 <- function(cdf, rim, extend=0.5)
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
.get_data3 <- function(cdf, rim, extend=0.5)
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
