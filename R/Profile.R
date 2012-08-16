Profile <- 
function(samples,Lib,peakData,r_thres=0.95, method = "dayNorm", minPairObs = 5){
    my.files <- RIfiles(samples)
    my.names <- sampleNames(samples)

    resInt <- Intensity(peakData)
    resRI  <- retIndex(peakData)
    resRT  <- retTime(peakData)
    libId  <- libId(Lib, sel = FALSE)

    # you can't set this paremeter lower than 5.
    minPairObs <- max(minPairObs, 5)

    # extract addictonal information
    addiData  <- libData(Lib)
    addiNames <- grep("^(Name|RI|Win_\\d*|SPECTRUM|TOP_MASS)$", colnames(addiData), value=TRUE, perl=TRUE, invert=TRUE)

    searchData <- data.frame(matrix(ncol=7,nrow=length(Lib)))
    rownames(searchData) <- 1:length(Lib)
    colnames(searchData) <- c("Mass_count", "Non_consecutive_Mass_count", "Sample_Count_per_Mass", "Masses",
        "RI", "Score_all_masses", "Score_cor_masses")

    medInt <- matrix(ncol=length(my.files),nrow=length(Lib))
    colnames(medInt) <- my.names
    rownames(medInt) <- 1:length(Lib)
    medRT <- medRI <- medInt

    res <- switch(method,
                dayNorm = dayNorm(samples, resInt),
                medianNorm = medianNorm(samples, resInt),
                none = resInt)
    res_log <- lapply(res,log2)

    options(warn=(-1))
    for(i in 1:length(Lib)){

        x <- as.character(selMass(Lib)[[i]])

 	    score_all <- NA
 	    score_cor <- NA

        # don't perform calculation with 1 mass
        if (length(x) > 1) {
 
            M <- res_log[[i]][x,]
  	        tmp <- cor(t(M), use="pair")
    	    # this counts the number of pair values that were used to calculate the correlation
	        # coefficient of every member of "tmp" and sets to 0 the cor.coef. with less than minPairObs
			tmp[is.finite(M) %*% is.finite(t(M)) <= minPairObs] <- 0

            # assume that the correlation of a metabolite with itself is always 1
            diag(tmp) <- 1

            tmp.max <- which.max(apply(tmp,1,function(x){ sum(x > r_thres, na.rm=T)}))
 
	        M <- res_log[[i]]
  	        tmp <- cor(t(M), use="pair")
 			tmp[is.finite(M) %*% is.finite(t(M)) < minPairObs] <- 0
            diag(tmp) <- 1
 
	        tmp.sel <- tmp[tmp.max,]
            y <- which(tmp.sel > r_thres)
  	  
 	  } else { # only 1 mass (no correlation)
			y <- 1 	  
 	  }

		# calculates scores if there are at least 3 top masses
		if(length(spectra(Lib)) > 0 & length(topMass(Lib)[[i]]) >= 3) {
			M <- apply(resInt[[i]], 1, median, na.rm = T)
			M[is.na(M)] <- 0
			score_all   <- Score(cbind(topMass(Lib)[[i]], M), spectra(Lib)[[i]])
			if(length(y) >= 3) {
				score_cor <- Score(cbind(topMass(Lib)[[i]][y], M[y]), spectra(Lib)[[i]])
			}
		}

    searchData$Mass_count[i]                 <- length(y)
    searchData$Non_consecutive_Mass_count[i] <- length(y) - sum(diff(sort(topMass(Lib)[[i]][y])) == 1)

        sampCountPerMass <- function(z) {
            tmp <- apply(z, 1, function(x) sum(is.finite(x)))
            tmp <- if(length(unique(tmp)) == 1) unique(tmp) else tmp
            paste(tmp, collapse=";")
        }

        searchData$Sample_Count_per_Mass[i] <- sampCountPerMass(res[[i]][y,,drop=FALSE])
    searchData$Masses[i] <- paste(topMass(Lib)[[i]][y], collapse=";")
    searchData$RI[i]     <- median(resRI[[i]][y,], na.rm=T)
    searchData$Score_all_masses[i] <- score_all
    searchData$Score_cor_masses[i] <- score_cor
 
		if(length(y) == 1) {
        medInt[i,] <- res[[i]][y,]
        medRI[i,] <- resRI[[i]][y,]
        medRT[i,] <- resRT[[i]][y,]
    } else {
        medInt[i,] <- apply(res[[i]][y,,drop=FALSE],2,median, na.rm=T)
        medRI[i,] <- apply(resRI[[i]][y,,drop=FALSE],2,median, na.rm=T)
        medRT[i,] <- apply(resRT[[i]][y,,drop=FALSE],2,median, na.rm=T)
    }
  }

    if(length(addiNames) > 0) {
        if(length(addiNames) == 1) {
            addiTemp <-  data.frame(addiData[,addiNames], stringsAsFactors=FALSE)
            colnames(addiTemp) <- addiNames
        } else {
            addiTemp <- addiData[,addiNames]
        }
        searchData <- cbind(searchData, addiTemp)
    }

  options(warn=0)
  return(new("tsProfile", peakData, info = cbind(Name = libName(Lib), Lib_RI = libRI(Lib), searchData),
    profInt=medInt, profRI=medRI, profRT=medRT))
}

# calculate scores. Inputs: two matrices with two columns. First and second
# columns are m/z and intensity, respectively. if "match" is TRUE, the score
# is obtain using the common masses. Otherwise, all the masses are used.

Score <- function(x, y, match = T) {
	if(all(x[,2] == 0) | all(y[,2] == 0))
		return(0)

	x <- x[order(x[,1]),]
	y <- y[order(y[,1]),]
	if(match) {
		x1 <- x[x[,1] %in% y[,1],2]
		y1 <- y[y[,1] %in% x[,1],2]
	} else {
		r <- range(x[,1], y[,1])
		x1 <- sapply(r[1]:r[2], function(z) if(any(x[,1] == z)) x[x[,1] == z,2] else 0)
		y1 <- sapply(r[1]:r[2], function(z) if(any(y[,1] == z)) y[y[,1] == z,2] else 0)
	}
	x1 <- 999 * x1 / max(x1)
	y1 <- 999 * y1 / max(y1)
	round((1 - sum(abs(x1-y1)) / sum(abs(x1+y1)))*1000)
}

