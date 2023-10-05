
`FindPeaks` <-
function(my.files, refLib, columns = NULL, showProgressBar = FALSE)
{
    my.names <- basename(my.files)
    if(is.list(refLib)) {
        # check that all matrices have the same size
        stopifnot(length(my.files) == length(refLib))
        stopifnot(all(sapply(refLib, ncol) == 3))
        nmz <- sapply(refLib,nrow)
        stopifnot(all(nmz == nmz[1]))
        nmz <- nmz[1]
        mz <- refLib[[1]][,2]
        if(is.null(rownames(refLib[[1]]))) {
            libId <- rep(1, nrow(refLib[[1]]))
        } else {
            libId <- rownames(refLib[[1]])
        }
    } else {
        nmz <- nrow(refLib)
        if(is.null(rownames(refLib))) {
            libId <- rep(1, nrow(refLib))
        } else {
            libId <- rownames(refLib)
        }
        mz    <- refLib[,2]
        refLib <- list(refLib)
    }

    resInt           <- matrix(nrow = nmz, ncol = length(my.files))
    colnames(resInt) <- my.names
    rownames(resInt) <- mz
    refLen           <- length(refLib)

    resRI <- resRT <- resInt

    if(showProgressBar)
        pb <- ProgressBar(title="Finding Peaks...", label="File in processing...")

    for (i in seq_along(my.files)) {
        lib <- refLib[ (i - 1) %% refLen + 1 ][[ 1 ]]
        out <- .c_find_peaks(
                            as.character(my.files[i]),   # MyFile
                            as.integer(mz),              # Mass
                            NULL,                        # RI_exp
                            as.numeric(lib[,1]),         # RI_min
                            as.numeric(lib[,3]),         # RI_max
                            FALSE,                       # useRT
                            'maxInt',                    # max intensity
                            columns)
        assert_that(!is.null(out),
                    msg=sprintf("Error processing file %s", my.files[i]))
        resInt[ out[[4]] + 1, i] <- out[[1]]
        resRI [ out[[4]] + 1, i] <- out[[2]]
        resRT [ out[[4]] + 1, i] <- out[[3]]

        if(showProgressBar)
            setProgressBar(pb, value=i/length(my.files),
                title=paste("Findind Peaks (", round(100*i/length(my.files)), "%)"),
                label=sprintf("Reading File %s", basename(my.files[i])))
    }

    if(showProgressBar)
        close(pb)

    # Transform resRI and resInt into a List
    mat2list <- function(x,g) lapply(unique(g), function(z) x[z==g,,drop=FALSE])
    resRI <- mat2list(resRI, libId)
    resRT <- mat2list(resRT, libId)
    resInt <- mat2list(resInt, libId)

    return(new("tsMSdata", RI = resRI, Intensity = resInt, RT = resRT))
}
