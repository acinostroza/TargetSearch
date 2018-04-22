# common functions

is_nullOrNA <- function(x)
{
    if(is.null(x)) return(TRUE)
    if(is.na(x)) return(TRUE)
    FALSE
}
