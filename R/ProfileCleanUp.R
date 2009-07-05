ProfileCleanUp <-
function(Profile, timeSplit=500, r_thres=0.95){

  if(nrow(profileInfo(Profile)) <= 1)
    return(Profile)

  o <- order(profileInfo(Profile)$RI)
  searchData2 <- profileInfo(Profile)[o,]
  medInt2 <-  Intensity(Profile)[o,]
  medRI2 <-  retIndex(Profile)[o,]
 
  if(sum(is.na(searchData2$RI) == FALSE) > 1) {
    medInt2 <-  medInt2[is.na(searchData2$RI) == F,]
    medRI2 <-  medRI2[is.na(searchData2$RI) == F,]
    searchData2 <- searchData2[is.na(searchData2$RI) == F,]
  } else if(sum(is.na(searchData2$RI) == FALSE) == 1) {
    medInt2 <- t( medInt2[is.na(searchData2$RI) == F,] )
    medRI2 <-  t( medRI2[is.na(searchData2$RI) == F,] )
    searchData2 <- searchData2[is.na(searchData2$RI) == F,]
    searchData2$final_sample_count <- sum(is.finite(medInt2))
    return(new("tsProfile", info = searchData2, Intensity = medInt2, RI = medRI2))    
  } else {
    return(Profile)
  }

  # as long as only consider MSTs with >=3 then quite robust against changing timeSplit
  # timeSplit <- 500
  
  tmGroups <- timeGroups(searchData2$RI, timeSplit)

  corGroups <- 1:nrow(searchData2)
  corData   <- searchData2[, colnames(searchData2) != "Lib_RI"]
  colnames(corData) <- paste("Cor_", colnames(corData), sep = "")
  corData$RI_dev <- format(searchData2$Lib_RI - searchData2$RI, digits = 3)
  corData$Cor_RI <- format(corData$Cor_RI, digits = 5)
  corData2 <- as.matrix(corData)

  for(i in 1:max(tmGroups)){
    if(sum(tmGroups == i) >1){
      tmp     <- cor( t(log2(medInt2[tmGroups == i,])), use = "pair")
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
  d <- duplicated(cbind(tmGroups,corGroups))
  MET_info <- cbind(searchData2[d == F,],corData[d == F,])
  MET <- medInt2[d == F,]
  MET_RI <- medRI2[d == F,]
  MET_info$final_sample_count <- apply(MET, 1, function(x) sum(is.finite(x)))
  rownames(MET) <- MET_info$Name
	rownames(MET_RI) <- MET_info$Name
  return(new("tsProfile", info = MET_info, Intensity = MET, RI = MET_RI))
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
