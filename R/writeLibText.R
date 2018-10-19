# function to write a tsLib object to a file.
# args: - lib: tsLib object
#       - file: a file

writeLibText <- function(lib, file)
{
	stopifnot(validObject(lib))
	libdata <- data.frame(Name=libName(lib), RI=libRI(lib),
				Win_1=RIdev(lib)[,1], Win_2=RIdev(lib)[,2],
				Win_3=RIdev(lib)[,3])
	libdata$SEL_MASS <- sapply(selMass(lib), function(x) paste(x, collapse=";"))
	libdata$TOP_MASS <- sapply(topMass(lib), function(x) paste(x, collapse=";"))
	if(length(spectra(lib)) > 0) {
		libdata$SPECTRUM <- sapply(spectra(lib), function(x)
					    paste(x[,1], x[,2], sep=":", collapse=" "))
	}
	extra <- !colnames(libData(lib)) %in% colnames(libdata)
	if(any(extra))
		libdata <- cbind(libdata, libData(lib)[, extra, drop=FALSE])
    write.table(libdata, file=file, row.names=FALSE, quote=FALSE, sep="\t")
}

