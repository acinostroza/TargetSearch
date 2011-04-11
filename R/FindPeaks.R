

FindPeaks <-
function(my.files, refLib, columns = c("SPECTRUM", "RETENTION_TIME_INDEX", "RETENTION_TIME"), showProgressBar = FALSE) {

    get.columns <- function(my.file, columns) {
        if(is.character(columns)) {
            header <- scan(my.file, what = "character", nlines = 1, quiet = TRUE)
            tmp <- sapply(columns, function(x) which( header == x ))
            if(length(unlist(tmp)) != length(columns))
                stop("Column name not found. Check your RI file.")
            columns <- unlist(tmp) - 1
        }
            return(columns)
    }

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
    }

    resInt           <- matrix(nrow = nmz, ncol = length(my.files))
    colnames(resInt) <- my.names
    rownames(resInt) <- mz

    resRI <- resRT <- resInt

    if(showProgressBar)
        pb <- ProgressBar(title="Finding Peaks...", label="File in processing...")
    for (i in 1:length(my.files)) {
        if(showProgressBar)
            setProgressBar(pb, value=i/length(my.files),
                title=paste("Findind Peaks (", round(100*i/length(my.files)), "%)"),
                label=sprintf("Reading File %s", basename(my.files[i])))

        cols <- get.columns(my.files[i], columns)
        out  <- if(is.list(refLib)) {
                    .Call("FindPeaks", as.character(my.files[i]), as.integer(refLib[[i]][,1]),
                        as.integer(refLib[[i]][,2]),	as.integer(refLib[[i]][,3]), as.integer(cols), PACKAGE="TargetSearch")
                } else { .Call("FindPeaks", as.character(my.files[i]), as.integer(refLib[,1]),
                        as.integer(refLib[,2]),	as.integer(refLib[,3]), as.integer(cols), PACKAGE="TargetSearch") }
        resInt[, i] <- out[[1]]
        resRI[, i]  <- out[[2]]
        resRT[, i]  <- out[[3]]
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
