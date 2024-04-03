#
# Functions to plot reference spectrum of an analyte/metabolite
#

#' Function to determine the m/z label positions
#'
#' this function is used to determine the positions of the labels in a m/z
#' spectrum, such that they not overlap with each other, nor they overlap with
#' the spectrum's vertical lines.
#'
#' @param x a numeric vector (the m/z values). These values should be integers;
#'      if not, they are internally converted to integers in order to compute
#'      their character width relative to the current plot.
#' @param y a numeric vector (the intensity values). Usually, the intensities
#'      are normalized from zero to thousand.
#' @param cex the magnification value for the labels
#' @param padding a one or two component vector representing a virtual "padding"
#'      for the labels in inches. The padding is required to determine whether
#'      two labels overlap with each others.
#' @param threshold a cutoff value: do not consider labels for absolute
#'      intensities smaller than this
#' @return A two component list. The first element is two-column matrix matrix
#'      representing the x, y coordinates of the labels, while the second are
#'      the respective labels.
#' @importFrom utils combn
#' @importFrom graphics strheight strwidth par
#' @importFrom assertthat assert_that is.number
`spectrum_label_pos` <- function(x, y, cex, padding, threshold)
{
    overlap <- function(a, b) {
        olap <- function(a, b, x, y) !(y <= a || x >= b)
        olap(a['x0'], a['x1'], b['x0'], b['x1']) &&
                                       olap(a['y0'], a['y1'], b['y0'], b['y1'])
    }

    assert_that(is.number(threshold))
    assert_that(length(padding) == 1 || length(padding) == 2)

    yo <- y ; xo <- x
    x <- x[ abs(y) > threshold]
    y <- y[ abs(y) > threshold]

    h <- strheight(x, cex=cex)
    w <- strwidth(x, cex=cex)
    u <- par('usr')
    p <- par('pin')
    sx <- (u[2] - u[1]) / p[1];
    sy <- (u[4] - u[3]) / p[2];
    padding <- padding * c(sx, sy)

    y1 <- y + sign(y) * (h + 2 * padding[2])
    rect <- cbind(x0= x - w/2 - padding[1], x1= x + w/2 + padding[1],
                  y0= pmin(y, y1), y1= pmax(y, y1), lab=floor(x))
    rect <- rect[order(-abs(y)),, drop=FALSE]
    vbar <- cbind(x0=xo, x1=xo, y0=pmin(0, yo), y1=pmax(0, yo))

    # remove rectangles that intersect with vertical bars
    ij <- expand.grid(i=seq(nrow(rect)), j=seq(nrow(vbar)))
    ov <- mapply(function(i, j) overlap(rect[i, ], vbar[j, ]), ij$i, ij$j)
    if(any(ov))
        rect <- rect[ -ij$i[ ov ],, drop=FALSE]

    # remove rectangles that intersect with each other
    if(nrow(rect) > 1) {
        ij <- data.frame(t(combn(nrow(rect), 2)))
        ov <- mapply(function(i, j) overlap(rect[i, ], rect[j, ]), ij$X1, ij$X2)
        if(any(ov))
            rect <- rect[ -ij$X2[ ov ],, drop=FALSE]
    }

    posx <- (rect[, 'x0'] + rect[, 'x1'])/2
    posy <- (rect[, 'y0'] + rect[, 'y1'])/2
    list(pos=cbind(posx, posy), labels=rect[, 'lab'])
}

#' Plot reference spectrum of a compound
#'
#' This function plots the reference spectram of an element of a `tsLib` object
#' see the manpage for more details.
#'
#' @param lib the library object to plot
#' @param libID the identifier of the metabolite to plot (defaults first one)
#' @param col the color of the vertical lines
#' @param main a main title for the plot
#' @param xlab a label for the x-axis, defaults to 'm/z'
#' @param ylab a label for the y-axis, defaults to 'Intensity'
#' @param label whether to display the m/z labels
#' @param cex the `cex` value used to scale the mz labels of the plot
#' @param label_col the color of the mz labels
#' @param sel_font the for of the selective masses
#' @param type the plot type. Only the value 'h' (vertical bars) is allowed.
#' @param ldots extra plotting parameters passed to plot
#' @return invisible
#' @importFrom assertthat assert_that is.count
`plotRefSpectra` <-
function(lib, libID, col='darkblue', main=NULL, xlab='m/z', ylab='Intensity',
         label=TRUE, cex=0.8, label_col=NULL, sel_font=4, type='h', ...)
{
    if(missing(libID))
        libID <- 1
    libID <- libID[1]
    assert_that(is.tsLib(lib))
    assert_that(is.flag(label))

    if(is.null(main))
        main <- libName(lib)[ libID ]

    if(is.null(spectra(lib)[[ libID ]]) || length(spectra(lib)[[ libID ]]) == 0) {
        warning('The spectrum is not available. Nothing to plot.')
        return(invisible())
    }
    x <- spectra(lib)[[ libID ]]
    m <- selMass(lib)[[ libID ]]

    assert_that(is.count(sel_font))

    plot(x, type='h', col=col, main=main, xlab=xlab, ylab=ylab, ...)

    pos <- spectrum_label_pos(x[,1], x[,2], cex, 0.04, 0.01 * max(x[,2]))
    if(label && !is.null(pos)) {
        labs <- pos$labels; pos <- pos$pos
        font <- rep(par('font'), nrow(pos))
        font[ pos[, 'posx'] %in% m ] <- sel_font
        text(pos, labels=labs, col=label_col, cex=cex, font=font)
    }
    return(invisible())
}

# vim: ts=4 sw=4 et:
