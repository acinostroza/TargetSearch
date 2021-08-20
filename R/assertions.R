
# assertions for checking object types

is.tsSample <- function(x) is(x, 'tsSample')
on_failure(is.tsSample) <- function(call, env) {
    paste0(deparse(call$x), " is not a `tsSample` object")
}

is.tsRim <- function(x) is(x, 'tsRim')
on_failure(is.tsRim) <- function(call, env) {
    paste0(deparse(call$x), " is not a `tsRim` object")
}

is.tsLib <- function(x) is(x, 'tsLib')
on_failure(is.tsLib) <- function(call, env) {
    paste0(deparse(call$x), " is not a `tsLib` object")
}

is.tsMSdata <- function(x) is(x, 'tsMSdata')
on_failure(is.tsMSdata) <- function(call, env) {
    paste0(deparse(call$x), " is not a `tsMSdata` object")
}

is.tsProfile <- function(x) is(x, 'tsProfile')
on_failure(is.tsProfile) <- function(call, env) {
    paste0(deparse(call$x), " is not a `tsProfile` object")
}

# useful assertions

is.sod <- function(x) length(x) == 1 || length(x) == 2
on_failure(is.sod) <- function(call, env) {
    paste0(deparse(call$x), " is neither a scalar nor its length is equal to 2")
}
