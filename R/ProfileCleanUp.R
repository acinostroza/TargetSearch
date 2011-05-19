ProfileCleanUp <-
function(Profile, timeSplit=500, r_thres=0.95, minPairObs=5){

  if(nrow(profileInfo(Profile)) <= 1)
    return(Profile)

  # you can't set this paremeter lower than 5.
  minPairObs <- max(minPairObs, 5)

  # It doesn't make sense to perform a clean up with few samples
  if(ncol(profileInt(Profile)) < minPairObs) {
    return(Profile)
  }

  o <- order(profileInfo(Profile)$RI)
  searchData2 <- profileInfo(Profile)[o,]
  medInt2     <- profileInt(Profile)[o,]
  medRI2      <- profileRI(Profile)[o,]
  medRT2      <- profileRT(Profile)[o,]
  intensity   <- Intensity(Profile)[o]
  RI          <- retIndex(Profile)[o]
  RT          <- retTime(Profile)[o]
 
  if(sum(is.na(searchData2$RI) == FALSE) > 1) {
    medInt2 <-  medInt2[is.na(searchData2$RI) == FALSE,]
    medRI2  <-  medRI2[is.na(searchData2$RI) == FALSE,]
    medRT2  <-  medRT2[is.na(searchData2$RI) == FALSE,]
    intensity   <- intensity[is.na(searchData2$RI) == FALSE]
    RI          <- RI[is.na(searchData2$RI) == FALSE]
    RT          <- RT[is.na(searchData2$RI) == FALSE]
    searchData2 <- searchData2[is.na(searchData2$RI) == FALSE,]

  } else if(sum(is.na(searchData2$RI) == FALSE) == 1) {
    medInt2 <- t( medInt2[is.na(searchData2$RI) == FALSE,] )
    medRI2 <-  t( medRI2[is.na(searchData2$RI) == FALSE,] )
    medRT2 <-  t( medRT2[is.na(searchData2$RI) == FALSE,] )
    intensity   <- intensity[is.na(searchData2$RI) == FALSE]
    RI          <- RI[is.na(searchData2$RI) == FALSE]
    RT          <- RT[is.na(searchData2$RI) == FALSE]
    searchData2 <- searchData2[is.na(searchData2$RI) == FALSE,]
    searchData2$final_sample_count <- sum(is.finite(medInt2))
    return(new("tsProfile", RI=RI, RT=RT, Intensity=intensity, info = searchData2,
     profInt = medInt2, profRI = medRI2, profRT = medRT2))
  } else {
    return(Profile)
  }

  # get additional columns
  foo <- which(colnames(searchData2) == "Score_cor_masses")
  addiInfo <- NULL
  if(length(foo) == 0) stop("Error in Profile. @info slot: column 'Score_cor_masses' doesn't exist")
  if(foo < ncol(searchData2)) {
    addiInfo <- colnames(searchData2)[(foo+1):ncol(searchData2)]
  }

  # as long as only consider MSTs with >=3 then quite robust against changing timeSplit
  # timeSplit <- 500
  
  tmGroups <- timeGroups(searchData2$RI, timeSplit)

  corGroups <- 1:nrow(searchData2)
  corData   <- searchData2[,! colnames(searchData2) %in% c("Lib_RI",addiInfo)]
  colnames(corData) <- paste("Cor_", colnames(corData), sep = "")
  corData$RI_dev <- format(searchData2$Lib_RI - searchData2$RI, digits = 3)
  corData$Cor_RI <- format(corData$Cor_RI, digits = 5)
  corData2 <- as.matrix(corData)

  for(i in 1:max(tmGroups)){
    if(sum(tmGroups == i) >1){
      M       <- log2(medInt2[tmGroups == i,])
      tmp     <- cor(t(M), use = "pair")
      tmp[is.finite(M) %*% is.finite(t(M)) < minPairObs] <- 0
      diag(tmp) <- 1
      tmp.max <- which.max(apply(tmp,1,function(x){ sum(x > r_thres, na.rm=T)}))
      tmp.sel <- tmp[tmp.max,]
      y <- which(tmp.sel > r_thres)
      if(length(y) > 1){
        corGroups[tmGroups == i][y] <- 1
        for(j in 1:ncol(corData))
	        corData2[tmGroups == i,j][y] <- paste(as.character(corData[tmGroups == i,j][y]), collapse=" | ")
      }
    }
  }
  o <- order(tmGroups,- corGroups, -searchData2[,3])
  tmGroups <- tmGroups[o]
  corGroups <- corGroups[o]
  searchData2 <- searchData2[o,]
  corData <- corData2[o,]
  medInt2 <- medInt2[o,]
  medRI2 <- medRI2[o,]
  medRT2 <- medRT2[o,]
  RI     <- RI[o]
  RT     <- RT[o]
  intensity <- intensity[o]
  
  d <- duplicated(cbind(tmGroups,corGroups))
  MET_info <- cbind(searchData2[d == F,],corData[d == F,])
  MET <- medInt2[d == F,]
  MET_RI <- medRI2[d == F,]
  MET_RT <- medRT2[d == F,]
  MET_info$final_sample_count <- apply(MET, 1, function(x) sum(is.finite(x)))
  rownames(MET) <- rownames(MET_RI) <- rownames(MET_RT) <- MET_info$Name
  RI <- RI[d == FALSE]
  RT <- RT[d == FALSE]
  intensity <- intensity[d == FALSE]
  return(new("tsProfile", RI=RI, RT=RT, Intensity=intensity, info = MET_info,
    profInt = MET, profRI = MET_RI, profRT = MET_RT))
}

timeGroups <- function(x, d) {
    if(length(x) == 1)
        return(1)
    if(length(x) == 0)
        return(NULL)
        
    o  <- order(x)
    gr <- numeric(length(o))
    for(i in 2:length(x))
        gr[i] <- if(is.na(x[o[i]]) | x[o[i]] - x[o[i-1]] > d) gr[i-1]+1 else gr[i-1]
    gr[order(o)]+1
}
