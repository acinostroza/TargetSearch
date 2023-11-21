# this function saves all results in a text format
`Write.Results` <-
function(Lib, metabProfile,
         quantMatrix=c('quantmass', 'maxint', 'maxobs', 'none'),
         prefix = NULL, selmass = FALSE)
{
    opt <- options(stringsAsFactors=FALSE)
    if(is.null_or_na(prefix))
        prefix <- paste("TargetSearch-", Sys.Date(), sep = "")
    assert_that(is.string(prefix))

    file.peak.int <- paste(prefix, ".peak.intensity.txt", sep ="")
    file.peak.ri  <- paste(prefix, ".peak.RI.txt", sep = "")
    file.pro.info <- paste(prefix, ".profile.info.txt", sep ="")
    file.pro.int  <- paste(prefix, ".profile.intensities.txt", sep ="")
    file.pro.ri   <- paste(prefix, ".profile.ri.txt", sep ="")
    file.quantmat <- paste(prefix, ".profile.quantmatrix.txt", sep="")

    qm <- match.arg(quantMatrix)
    if(qm != 'none') {
        M <- quantMatrix(Lib, metabProfile, value=qm, selmass=selmass)
        M <- data.frame(libID=rownames(M), Name=attr(M, "libNames"),
                        quantMass=attr(M, "quantMass"), isSelectiveMass=attr(M, "isSelMass"),
                        isCorMass=attr(M, "isCorMass"), M)
        write.table(M, file=file.quantmat, sep="\t", quote=FALSE, row.names=FALSE)
    }

    write.table(profileInfo(metabProfile), file = file.pro.info, sep = "\t", quote = FALSE)
    write.table(profileInt(metabProfile), file = file.pro.int, sep = "\t", quote = FALSE)
    write.table(profileRI(metabProfile), file = file.pro.ri, sep = "\t", quote = FALSE)

    metInfo   <- profileInfo(metabProfile)
    metPrf    <- as.list(metabProfile)
    libIndex  <- attr(metPrf$RI, "index")
    libId     <- rownames(metInfo)
    med_RI    <- format(apply(metPrf$RI, 1, median, na.rm = TRUE), digits = 3)

    Name      <- as.character(metInfo[libIndex, "Name"])
    libRI     <- metInfo[libIndex, "Lib_RI"]
    mass      <- unlist(topMass(Lib)[libId])
    is_sel    <- unlist(lapply(libId, function(x) topMass(Lib)[[x]] %in% selMass(Lib)[[x]]))
    corMass   <- lapply(strsplit(metInfo$Masses, ";"), function(x) as.numeric(x))
    is_cor    <- unlist(lapply(seq(libId), function(x) topMass(Lib)[[libId[x]]] %in% corMass[[x]]))
    
    Out <- data.frame(libId=libIndex, Name= Name, Lib_RI=libRI, Mass = mass,
         med_RI, IS_SEL = is_sel, IS_COR=is_cor)

    write.table( data.frame(Out, metPrf$Intensity, row.names = NULL, check.names = FALSE),
                file = file.peak.int, sep = "\t", quote = FALSE, row.names = FALSE)
    write.table( data.frame(Out, metPrf$RI, row.names = NULL, check.names = FALSE),
                file = file.peak.ri, sep = "\t", quote = FALSE, row.names = FALSE)
    options(opt)
    invisible()
}

# function to create a MSP file that can be viewed with NIST
# metlib : library
# metprof: at tsProfile object (output of Profile or ProfileCleanUp)
# file: msp output file name
# append: should data be appended to 'file'

writeMSP <- function (metlib, metprof, file, append = FALSE)
{
    if (is.logical(append) == FALSE) {
        stop("append must be logical")
    }
    tmp <- Intensity(metprof)
    # check ids
    if(!all(names(tmp) %in% rownames(libData(metlib))))
        stop("library identifiers don't match profile identifiers.")

    file <- file(file, ifelse(append[1], "a", "w"))
    for (id in names(tmp)) {
        x <- tmp[[id]]
        i <- which(id == rownames(libData(metlib)))
        if (all(is.na(x))) {
            next
        }
        x.median <- apply(x, 1, median, na.rm = TRUE)
        x.median[is.na(x.median)] <- 0
        x.median <- round(999 * x.median/max(x.median))
        mz <- topMass(metlib)[[id]][x.median > 0]
        x.median <- x.median[x.median > 0]
        if(length(x.median) < 2)
            next
        cat(sprintf("Name: %s", libName(metlib)[i]), file = file,
            sep = "\n")
        cat(sprintf("Synon: RI: %.1f", medRI(metlib)[i]), file = file,
            sep = "\n")
        cat(sprintf("Synon: MST SEL MASS: %s", substring(paste(selMass(metlib)[[id]],
            collapse = ";"), 1, 230)), file = file, sep = "\n")
        cat(sprintf("Num Peaks: %d", length(x.median)), file = file,
            sep = "\n")
        foo <- paste(mz, " ", x.median, ";", sep = "")
        foo <- split(foo, rep(1:ceiling(length(foo)/6), each = 6)[1:length(foo)])
        cat(sapply(foo, function(x) paste(x, collapse = " ")),
            file = file, sep = "\n")
        cat("", file = file, sep = "\n")
    }
    close(file)
    invisible()
}

# vim: set ts=4 sw=4 expandtab:
