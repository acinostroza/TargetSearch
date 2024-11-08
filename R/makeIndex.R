#' Make an index vector out of selective/top masses
#'
#' @param lib a tsLib object.
#' @param sel a logical flag. If `TRUE` return an integer index for selective
#'        masses, otherwise for top masses.
#' @return an integer vector representing indices
#' @note This replaces the method `libId`
#'
`makeIndex` <- function(lib, sel=TRUE)
{
    assert_that(is.tsLib(lib), is.flag(sel))
    mass <- if(sel) selMass(lib) else topMass(lib)
    rep.int(seq(lib), vapply(mass, length, 0))
}

# vim: set ts=4 sw=4 et:
