# common functions

is_nullOrNA <- function(x)
{
    if(is.null(x)) return(TRUE)
    if(any(is.na(x))) return(TRUE)
    FALSE
}

#' like `dirname`, but checks if all dirs are the same
#'
#' @param x a vector of file paths
#' @return the dirname for each member of x. If all names are equal, return
#'         the unique value
.dirname <- function(x)
{
    ret <- dirname(x)
    unq <- unique(ret)
    if(length(unq) == 1)
        return(unq)
    return(ret)
}

#' update file path
#'
#' Updates the path of new files. If the length of new is 1, then apply the
#' the same path to all files. if NULL or NA, return old paths
#' Assumes that both length are equal or the paths can be recycled. Raises
#' and error if it is not possible
#'
#' @param old a vector of file paths
#' @param new a vector of new paths
#' @return the updated file paths
.setpath <- function(old, new)
{
    if(is_nullOrNA(new))
        return(old)
    if(length(old) %% length(new) != 0)
        stop("lengths must be multiples of each other")
    if(length(new) > length(old)) {
        warning('length of `new` is longer than `old`. Extra elements will be ignored')
        new <- new[seq(length(old))]
    }
    new <- rep(new, length(old)/length(new))
    file.path(new, basename(old))
}

#' trim file extension
#'
#' remove known file extensions from file path. if none found, then does
#' not do anything.
#'
#' @param files file names
#' @param exts  file extensions to trim. the extensions are converted to
#'              to regular expressions
#' @return files with extension trimmed.
.trim_file_ext <- function(files, exts)
{
    e <- paste(exts, collapse="|")
    e <- sprintf("\\.(%s)$", e)
    str_remove(files, regex(e, ignore_case=TRUE))
}

#' add file extension
#'
#' add file extensions to extension trimmed file names and checks if
#' the file exists. In case of ambiguity, takes the first
#' @param a list of file paths
#' @param a list of file extensions
#' @return take the first path/extension combo that exists, if none exists
#'         take the first.
.add_file_ext <- function(files, exts)
{
    f <- lapply(exts, function(e) paste(files, e, sep="."))
    f <- lapply(seq(files), function(i) sapply(f, getElement, i))
    sapply(f, function(z) {
               y <- file.exists(z)
               if(all(!y)) z[1] else z[y][1] })
}

# vim: set ts=4 sw=4 et:
