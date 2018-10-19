# low level search of peaks in a RI file.
`getAllPeaks` <-
function(file, ref, useRT=FALSE, searchType = c("all", "minRI", "maxInt"),
         columns = c("SPECTRUM", "RETENTION_TIME_INDEX", "RETENTION_TIME"))
{
    if(!all(c('mz', 'minRI', 'maxRI') %in% colnames(ref)))
        stop("Error: missing columns in 'ref'")

    opt <- as.integer(get.file.format.opt(file, columns))
    searchType <- pmatch(searchType, c("all", "minRI", "maxInt"))
    z <- .Call(c_find_peaks, file, as.integer(ref[,'mz']), NULL,
               as.numeric(ref[,'minRI']), as.numeric(ref[,'maxRI']), opt, useRT, searchType)
    z <- do.call('cbind', z)
    z[,4] <- z[,4] + 1
    colnames(z) <- c('Int', 'RI', 'RT', 'rowid')
    cbind(z, mz=ref[z[, 'rowid'],2])
}

`FindAllPeaks` <-
function(samples, Lib, libID, dev=NULL, mz=NULL, RI=NULL,
         mz_type = c('selMass', 'quantMass', 'topMass'),
         columns = c("SPECTRUM", "RETENTION_TIME_INDEX", "RETENTION_TIME"))
{
    if(is_nullOrNA(RI)) #
        RI <- if(!is.na(medRI(Lib)[libID])) medRI(Lib)[libID] else libRI(Lib)[libID]

    if(is_nullOrNA(dev))
        dev <- RIdev(Lib)[libID, 1]

    if(is_nullOrNA(mz)) {
        method <- switch(match.arg(mz_type),
                         selMass=selMass, quantMass=quantMass, topMass=topMass)
        mz <- method(Lib)[[libID]]
    }

    ref   <- cbind(minRI=RI-dev, mz=mz, maxRI=RI + dev)
    peaks <- lapply(RIfiles(samples), getAllPeaks, ref, columns=columns)

    n  <- rep.int(1:length(peaks), sapply(peaks, nrow))
    pk <- do.call('rbind', peaks)
    pk <- cbind(pk, fid=n)

    if(nrow(pk) == 0)
        return(NULL)

    rownames(pk) <- NULL
    return(pk[, c('Int', 'RI', 'RT', 'mz', 'fid'), drop=FALSE])
}

# vim: set ts=4 sw=4 et:
