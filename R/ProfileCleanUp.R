ProfileCleanUp <-
function(Profile, timeSplit=500, r_thres=0.95){
  o <- order(profileInfo(Profile)$RI)
  searchData2 <- profileInfo(Profile)[o,]
  medInt2 <-  Intensity(Profile)[o,]
  medInt2 <-  medInt2[is.na(searchData2$RI) == F,]
  medRI2 <-  retIndex(Profile)[o,]
  medRI2 <-  medRI2[is.na(searchData2$RI) == F,]
  searchData2 <- searchData2[is.na(searchData2$RI) == F,]
  all_cor2 <- cor(t(log2(medInt2)), use="pair")
                               	
  # as long as only consider MSTs with >=3 then quite robust against changing timeSplit
  # timeSplit <- 500
  
  timeGroups <- numeric(length=nrow(searchData2))
  tmp <- diff(searchData2$RI)
  for(i in 2:nrow(searchData2)){
    if(tmp[i-1] <= timeSplit) {
      timeGroups[i] <- timeGroups[i-1]
    }
    else {
      timeGroups[i] <- timeGroups[i-1] + 1
    }
  }
  timeGroups <- timeGroups + 1

  corGroups <- 1:nrow(searchData2)
  corData   <- searchData2[, colnames(searchData2) != "Lib_RI"]
  colnames(corData) <- paste("Cor_", colnames(corData), sep = "")
  corData$RI_dev <- searchData2$Lib_RI - searchData2$RI
  corData2 <- as.matrix(corData)

  for(i in 1:max(timeGroups)){
    if(sum(timeGroups == i) >1){
      tmp <- all_cor2[timeGroups == i, timeGroups == i]
      tmp.max <- which.max(apply(tmp,1,function(x){ sum(x > r_thres, na.rm=T)}))
      tmp.sel <- tmp[tmp.max,]
      y <- which(tmp.sel > r_thres)
      if(length(y) > 1){
        corGroups[timeGroups == i][y] <- 1
        for(j in 1:ncol(corData))
	        corData2[timeGroups == i,j][y] <- paste(as.character(corData[timeGroups == i,j][y]), collapse=" | ")
      }
    }
  }
  o <- order(timeGroups,- corGroups, -searchData2[,3])
  timeGroups <- timeGroups[o]
  corGroups <- corGroups[o]
  searchData2 <- searchData2[o,]
  corData <- corData2[o,]
  medInt2 <- medInt2[o,]
  medRI2 <- medRI2[o,]
  d <- duplicated(cbind(timeGroups,corGroups))
  MET_info <- cbind(searchData2[d == F,],corData[d == F,])
  MET <- medInt2[d == F,]
  MET_RI <- medRI2[d == F,]
  MET_info$final_sample_count <- apply(MET, 1, function(x) sum(is.finite(x)))
  rownames(MET) <- MET_info$Name   
	rownames(MET_RI) <- MET_info$Name
  return(new("tsProfile", info = MET_info, Intensity = MET, RI = MET_RI))
}

