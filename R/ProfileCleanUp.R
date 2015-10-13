`ProfileCleanUp` <-
function(Profile, timeSplit=500, r_thres=0.95, minPairObs=5,
    prioritization=c('mass','score'), corMass=1, score=0,
    show=c('unidentified','knowns','full'))
{
  if(nrow(profileInfo(Profile)) <= 1)
    return(Profile)

  # you can't set this paremeter lower than 5.
  if(minPairObs < 5)
    warning("'minPairObs' cannot be set to a value lower than 5. ",
      "Using 5 instead.")
  minPairObs <- max(minPairObs, 5)

  # It doesn't make sense to perform a clean up with few samples
  if(ncol(profileInt(Profile)) < minPairObs) {
    warning("The number of samples is less than 'minPairObs'")
    return(Profile)
  }

  opt <- options(stringsAsFactors=FALSE)

  prioritization <- match.arg(prioritization)
  show           <- match.arg(show)

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
    options(opt)
    return(new("tsProfile", RI=RI, RT=RT, Intensity=intensity, info = searchData2,
     profInt = medInt2, profRI = medRI2, profRT = medRT2))
  } else {
    options(opt)
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

  corGroups <- rep.int(1, nrow(searchData2))
  corData   <- searchData2[,! colnames(searchData2) %in% c("Lib_RI",addiInfo)]
  colnames(corData) <- paste("Cor_", colnames(corData), sep = "")
  corData$RI_dev <- format(searchData2$Lib_RI - searchData2$RI, digits = 3)
  corData$Cor_RI <- format(corData$Cor_RI, digits = 5)

  opt_warn <- options(warn=-1) # suppress correlation warnings
  for(i in 1:max(tmGroups)) {
    if(sum(tmGroups == i) > 1) {
      M       <- log2(medInt2[tmGroups == i,])
      tmp     <- cor(t(M), use = "pair")
      tmp[is.finite(M) %*% is.finite(t(M)) < minPairObs] <- 0
      tmp[is.na(tmp) | is.nan(tmp) ] <- 0
      diag(tmp) <- 1
      hc <- hclust(as.dist(1-tmp), 'single')
      ct <- cutree(hc, h=1-r_thres)
      corGroups[tmGroups == i] <- ct
    }
  }
  options(opt_warn)

  Groups <- cbind(Time_group=tmGroups, Cor_group=corGroups,
               searchData2[,c("Mass_count","Score_all_masses")])
  o <- rankGroups(Groups, prioritization, corMass, score)
  Groups   <- Groups[o,]
  tmGroups <- tmGroups[o]
  corGroups <- corGroups[o]
  searchData2 <- searchData2[o,]
  corData <- corData[o,]
  medInt2 <- medInt2[o,]
  medRI2 <- medRI2[o,]
  medRT2 <- medRT2[o,]
  RI     <- RI[o]
  RT     <- RT[o]
  intensity <- intensity[o]

  # create the cor_data columns according to the prioritization order
  corData2 <- as.matrix(corData)
  for(i in 1:max(tmGroups)) {
    if(sum(tmGroups == i) > 1) {
      ct <- corGroups[tmGroups == i]
      for(ii in 1:max(ct)) {
        if(sum(ct == ii) > 1) {
          for(j in 1:ncol(corData))
	          corData2[tmGroups == i,j][ct == ii] <-
              paste(as.character(corData[tmGroups == i,j][ct == ii]), collapse=" | ")
        }
      }
    }
  }
  corData <- corData2

  # identify the bad candidates
  k <- which(searchData2$Mass_count < corMass | searchData2$Score_all_masses < score)
  if(!is.character(searchData2$Name))
    searchData2$Name <- as.character(searchData2$Name)
  searchData2$Name[k] <- sprintf("Unidentified %s [%s]",
                            format(searchData2$RI[k], digits=6), searchData2$Name[k])

  d <- duplicated(cbind(tmGroups,corGroups))
  if(show == 'knowns') {
    d[k] <- TRUE
  } else if(show == 'full') {
    d[]  <- FALSE
    searchData2 <- cbind(searchData2, Groups[, 1:2])
  }
  MET_info <- cbind(searchData2[d == F,],corData[d == F,])
  MET <- medInt2[d == F,]
  MET_RI <- medRI2[d == F,]
  MET_RT <- medRT2[d == F,]
  MET_info$final_sample_count <- apply(MET, 1, function(x) sum(is.finite(x)))
  rownames(MET) <- rownames(MET_RI) <- rownames(MET_RT) <- MET_info$Name
  RI <- RI[d == FALSE]
  RT <- RT[d == FALSE]
  intensity <- intensity[d == FALSE]

  options(opt)
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

# order the metabolites according to the prioritization parameter
# score and correlating masses considering the thresholds
rankGroups <- function(x, p, cm, sc)
{
    stopifnot(ncol(x) == 4)
    x <- cbind(x, x[,3] >= cm & x[,4] >= sc)
    if(p == 'mass')
        o <- order(x[,1], x[,2], -x[,5], -x[,3], -x[, 4])
    else
        o <- order(x[,1], x[,2], -x[,5], -x[,4], -x[, 3])
}

# vim: set ts=4 sw=4 expandtab:
