#' plot peak from RI files
#'
#' @param obj list of files of tsSample object
#' @param massValues the m/z values (two-element vector)
#' @param timeRange the time range (two-element vector)
#' @param useRT use RT as time units (RI if FALSE)
#' @param showRT use RT as y-axis (RI if FALSE)
#' @param sizefun funtion to scale the plot size
#' @param plot should a plot be drawn?
#' @param ... extra plotting parameters
#' @return a matrix with the data used for plotting
`ri_plot_peak` <-
function(obj, massValues, timeRange, useRT=TRUE, showRT=useRT, sizefun=NULL, plot=TRUE, ...)
{
    if(inherits(obj, 'tsSample')) {
        ri <- RIfiles(obj)
    } else if(is.character(obj)) {
        ri <- obj
    } else {
        stop('Object must inherit `tsSample` or be a character vector')
    }
    assert_that(is.null(sizefun) || is.function(sizefun))
    assert_that(is.flag(useRT))
    assert_that(is.flag(showRT))
    assert_that(is.flag(plot))

    size_fun <- function(z) pmax(0.1, (log10(z) - 2) * 9 / 8 + 1 / 2)

    # i <- =mz
    # j <- =sample
    optfun <- function(i, j, ...) {
        S <- function(x) paste0('s', x)
        J <- function(x, i) x[ ((i - 1) %% length(x)) + 1]
        opts <- list(...)
        optnames <- c('col', 'pch', 'bg', 'cex')
        for(name in optnames) {
            if(!is.null(opts[[ S(name) ]])) {
                opts[[ name ]] <- J(opts[[ S(name) ]], j)
                opts[[ S(name) ]] <- NULL
            } else if(!is.null(opts[[ name ]])) {
                opts[[ name ]] <- J(opts[[ name ]], i)
            }
        }
        opts[!duplicated(names(opts))]
    }

    dat <- lapply(ri, ri_data_extract, massValues, timeRange, useRT)
    n <- sapply(dat, nrow)
    dat <- cbind(do.call(rbind, dat), sample=rep(seq_along(dat), n))
    if(nrow(dat) == 0) {
        warning("Nothing to plot")
        return(invisible())
    }
    rownames(dat) <- NULL

    if(!plot)
        return(invisible(dat))

    time <- if(showRT) 'RT' else 'RI'
    mzid <- match(dat[,'mz'], massValues)
    opts <- optfun(mzid, dat[, 'sample'], x=dat[, c('sample', time), drop=FALSE], ..., type='p')

    if(is.null(opts$cex)) {
        fun <- if(is.null(sizefun)) size_fun else sizefun
        opts$cex <- fun(dat[, 'Intensity'])
    }

    do.call('plot', opts)

    return(invisible(dat))
}
