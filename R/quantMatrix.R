quantMatrix <- function(Lib, metabProfile, value = "maxint") {
    value <- pmatch(value, c("maxint", "maxobs"))
    Int   <- Intensity(metabProfile)
    qM    <- quantMass(Lib)
    if(length(qM) == 0)
        qM <- numeric(length(Lib))
    M <- matrix(nrow=length(Int),ncol=ncol(Int[[1]]),
        dimnames = list(profileInfo(metabProfile)$Name,colnames(profileInt(metabProfile))))

    id <- rownames(profileInfo(metabProfile))
    stopifnot(all(id %in% names(selMass(Lib))))
    id <- match(id, names(selMass(Lib)))

    sM <- profileInfo(metabProfile)$Masses
    attr(M, "quantMass") <- character(length(id))
    attr(M, "isSelMass") <- logical(length(id))
    attr(M, "isCorMass") <- logical(length(id))

    for(i in 1:length(id)) {
        if(is.na(qM[id[i]]) | qM[id[i]] == 0) {
            if(is.na(sM[i]) | sM[i] == "") {
                attr(M, "isSelMass")[i] <- FALSE
                attr(M, "quantMass")[i] <- NA
                attr(M, "isCorMass")[i] <- FALSE
                next
            }
            mz <- unlist(strsplit(sM[i], ";"))
            attr(M, "isCorMass")[i] <- TRUE
            int <- Int[[i]][mz,]
            mz.sel <- mz %in% selMass(Lib)[[id[i]]]
            if(any(mz.sel)) {
                if(!is.null(dim(int))) int <- int[mz.sel,]
                mz  <- mz[mz.sel]
                attr(M, "isSelMass")[i] <- TRUE
            } else {
                attr(M, "isSelMass")[i] <- FALSE
            }
            if(is.null(dim(int))) {
                M[i,] <- int
                attr(M, "quantMass")[i] <- mz
                next
            }
            mi <- which.max(apply(int,1,median,na.rm=TRUE))
            mc <- which.max(apply(int,1,function(x) sum(is.na(x)==FALSE)))
            if(value == 1) {
                M[i,] <- int[mi,]
                attr(M, "quantMass")[i] <- mz[mi]
            } else if(value == 2) {
                M[i,] <- int[mc,]
                attr(M, "quantMass")[i] <- mz[mc]
            }
        } else {
            if(any(qM[id[i]] == selMass(Lib)[[id[i]]])) {
                M[i,] <- Int[[i]][as.character(qM[id[i]]),]
                attr(M, "quantMass")[i] <- qM[id[i]]
                attr(M, "isSelMass")[i] <- TRUE
            } else {
                attr(M, "quantMass")[i] <- qM[id[i]]
                attr(M, "isSelMass")[i] <- FALSE
                message('Warning: mz=', qM[id[i]], ' not found in metabolite \"', libName(Lib)[id[i]], '"')
            }
            mz <- unlist(strsplit(sM[i], ";"))
            attr(M, "isCorMass")[i] <- any(qM[id[i]] %in% mz)
        }
    }
    M
}
