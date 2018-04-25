# common functions

is_nullOrNA <- function(x)
{
    if(is.null(x)) return(TRUE)
    if(any(is.na(x))) return(TRUE)
    FALSE
}
