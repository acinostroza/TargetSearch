
`FindPeaks` <-
function(my.files, refLib, columns = c("SPECTRUM", "RETENTION_TIME_INDEX", "RETENTION_TIME"), showProgressBar = FALSE)
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
    }

    resInt           <- matrix(nrow = nmz, ncol = length(my.files))
    colnames(resInt) <- my.names
    rownames(resInt) <- mz

    resRI <- resRT <- resInt

    if(showProgressBar)
        pb <- ProgressBar(title="Finding Peaks...", label="File in processing...")

    for (i in 1:length(my.files)) {

        # Guess the file format. See file "file.R"
        opts <- get.file.format.opt(my.files[i], columns)
        # set option too select for max intensity (default behaviour)
        searchType <- pmatch("maxInt", c("all", "minRI", "maxInt"))

        out  <- if(is.list(refLib)) {
                    .Call(c_find_peaks,
                            as.character(my.files[i]),   # MyFile
                            as.integer(mz),              # Mass
                            NULL,                        # RI_exp
                            as.numeric(refLib[[i]][,1]), # RI_min
                            as.numeric(refLib[[i]][,3]), # RI_max
                            as.integer(opts),            # Options
                            FALSE,                       # useRT
                            searchType,                  # max intensity
                            PACKAGE="TargetSearch")
                } else {
                    .Call(c_find_peaks,
                            as.character(my.files[i]), # MyFile
                            as.integer(mz),            # Mass
                            NULL,                      # RI_exp
                            as.numeric(refLib[,1]),    # RI_min
                            as.numeric(refLib[,3]),    # RI_max
                            as.integer(opts),          # Options
                            FALSE,                     # useRT
                            searchType,                # max intensity
                            PACKAGE="TargetSearch")
                }
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
