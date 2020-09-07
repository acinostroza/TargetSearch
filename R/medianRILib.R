# ChangeLog:
# 07.10.2008: the function was change to support the new Library format.

medianRILib <- function(samples, Lib, makeReport = FALSE, pdfFile = "medianLibRep.pdf",
	columns = NULL, showProgressBar = FALSE) {

	my.files   <- RIfiles(samples)
	refLib     <- refLib(Lib, w = 1, sel = TRUE)
	libId      <- libId(Lib, sel = TRUE)
	resPeaks   <- FindPeaks(my.files, refLib, columns, showProgressBar)
	med_RI     <- sapply(retIndex(resPeaks), median, na.rm = T)
	medRI(Lib) <- med_RI

	if(makeReport == TRUE)
	 	plotAllRIdev(Lib, resPeaks, pdfFile)

  return(Lib)
}

