Profile <- 
function(samples,Lib,peakData,r_thres=0.95, method = "dayNorm", minPairObs = 5){
	my.files <- RIfiles(samples)
	my.names <- sampleNames(samples)

  resInt <- Intensity(peakData)
  resRI  <- retIndex(peakData)
  libId  <- libId(Lib, sel = FALSE)
  
  # extract addictonal information
  addiData  <- libData(Lib)
  addiNames <- grep("^(Name|RI|Win_\\d*|SPECTRUM|TOP_MASS)$", colnames(addiData), value=TRUE, perl=TRUE, invert=TRUE)

  searchData <- data.frame(matrix(ncol=7,nrow=length(Lib)))
  rownames(searchData) <- 1:length(Lib)
  colnames(searchData) <- c("Mass_count", "Non_consecutive_Mass_count", "Sample_count", "Masses",
		"RI", "Score_all_masses", "Score_cor_masses")
		
  medInt <- matrix(ncol=length(my.files),nrow=length(Lib))
  colnames(medInt) <- my.names
  rownames(medInt) <- 1:length(Lib)
  medRI <- medInt
  
  res <- switch(method, 
						dayNorm = dayNorm(samples, resInt),
						medianNorm = medianNorm(samples, resInt),
						none = resInt)
  res_log <- log2(res)
  
  options(warn=(-1))
  for(i in 1:length(Lib)){

    x     <- which(libId == i)
    x.sel <- x[1:length(selMass(Lib)[[i]])]
    
 	  score_all <- NA
 	  score_cor <- NA

    # don't perform calculation with 1 mass
    if (length(x) > 1) {
 
	    M <- res_log[x.sel,]
  	  tmp <- cor(t(M), use="pair")
    	# this counts the number of pair values that were used to calculate the correlation
	    # coefficient of every member of "tmp" and sets to 0 the cor.coef. with less than minPairObs 
			tmp[is.finite(M) %*% is.finite(t(M)) <= minPairObs] <- 0
    
    	tmp.max <- which.max(apply(tmp,1,function(x){ sum(x > r_thres, na.rm=T)}))
    
	    M <- res_log[x,]
  	  tmp <- cor(t(M), use="pair")
 			tmp[is.finite(M) %*% is.finite(t(M)) < minPairObs] <- 0
    
	    tmp.sel <- tmp[tmp.max,]
    
  	  y <- which(tmp.sel > r_thres)
  	  
 	  } else { # only 1 mass (no correlation)
			y <- 1 	  
 	  }

		# calculates scores 	  
 	  if(length(spectra(Lib)) > 0 & length(x) >= 3) {
 	  	M <- apply(resInt[x, ], 1, median, na.rm = T)
 	  	M[is.na(M)] <- 0
 	  	score_all   <- Score(cbind(topMass(Lib)[[i]], M), spectra(Lib)[[i]])
 	  	if(length(y) >= 3) {
 	  		score_cor <- Score(cbind(topMass(Lib)[[i]][y], M[y]), spectra(Lib)[[i]])
 	  	}
 	  }
 	  
 	  
    searchData$Mass_count[i]                 <- length(y)
    searchData$Non_consecutive_Mass_count[i] <- length(y) - sum(diff(sort(topMass(Lib)[[i]][y])) == 1)
		searchData$Sample_count[i]               <- sum(is.finite(res[x[y],]))
    searchData$Masses[i] <- paste(topMass(Lib)[[i]][y], collapse=";")
    searchData$RI[i]     <- median(resRI[x[y],], na.rm=T)
    searchData$Score_all_masses[i] <- score_all
    searchData$Score_cor_masses[i] <- score_cor
    
		if(length(y) == 1) {
        medInt[i,] <- res[x[y],]
        medRI[i,] <- resRI[x[y],]
    } else {
        medInt[i,] <- apply(res[x[y],],2,median, na.rm=T)
        medRI[i,] <- apply(resRI[x[y],],2,median, na.rm=T)
    }
  }

  if(length(addiNames) > 0) searchData <- cbind(searchData, addiData[,addiNames])

  options(warn=0)
	tmp <- new("tsMSdata", RI = medRI, Intensity = medInt)  
  return(new("tsProfile", tmp, info = cbind(Name = libName(Lib), Lib_RI = libRI(Lib), searchData)))
}

# calculate scores. Inputs: two matrices with two columns. First and second
# columns are m/z and intensity, respectively. if "match" is TRUE, the score
# is obtain using the common masses. Otherwise, all the masses are used.

Score <- function(x, y, match = T) {
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
