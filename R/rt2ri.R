`rt2ri` <-
function(rtTime, observed, standard) {

  ##<..Beg Rdocu..>
  ## ~name~
  ##   rt2ri
  ## ~title~
  ##   RT to RI
  ## ~description~
  ##   Convert retention times to retention indices based on observed
  ##   FAME RI and their standard values
  ## ~usage~
  ##   rt2ri(rtTime, observed, standard=c(262320, 323120, 381020,
  ##   487220, 582620, 668720, 747420, 819620, 886620, 948820,
  ##   1006900, 1061700, 1113100))
  ## ~arguments~
  ##   ~-rtTime~
  ##     The extracted RT's to convert
  ##   ~-observed~
  ##     The observed FAME RT's
  ##   ~-standard~
  ##     The standard RI for each FAME
  ## ~details~
  ##   Linear interpolation, interpolation outside bounds are done
  ##   like Jan's scripts with continued linear interpolation from the
  ##   last two FAME's
  ## ~value~
  ##   The converted RI
  ## ~seealso~
  ##   AttachRI from Jan's scripts
  ## ~examples~
  ##   bla <- rt2ri(1:100, rnorm(13))
  ## ~author~
  ##   Henning Redestig <redestig[at]mpimp-golm.mpg.de>
  ##>..End Rdocu..<

  if(length(observed) != length(standard))
   	stop("Function rt2ri: 'observed' and 'standard' parameters have different lengths")

  lost <- is.na(observed)
  if(any(lost)) {
    message("Missing RI-markers")
    observed <- observed[!lost]
    standard <- standard[!lost]
  }
	os <- order(standard)
	standard <- standard[os]
	observed <- observed[os]

  retData <- data.frame(obs=observed, sta=standard)
  o <- order(rtTime)
  revo <- match(1:length(o), o)

  ## approximate
  appr <- approx(observed, standard, xout=rtTime[o])
  ## fix the outside boundaries NA's
  belowLm <- lm(sta~obs, retData, subset=1:2)
  aboveLm <- lm(sta~obs, retData, subset=(length(standard) - 1):length(standard))
  ## get the indices of the retention times that could not be approximated
  below <- which(is.na(appr$y[1:(length(appr$y) / 2)]))
  above <- which(is.na(appr$y[(floor(length(appr$y) / 2)):length(appr$y)])) + floor(length(appr$y) / 2) - 1
  appr$y[below] <- predict(belowLm, new=list(obs=appr$x[below]))
  appr$y[above] <- predict(aboveLm, new=list(obs=appr$x[above]))
  appr$y[revo]
}

#########################################################################################

# Function to transform a matrix of retention indices to a matrix of retention time
# usage: ri2rt(x, rt.observed, ri.standard)
# args:
#      x:           A matrix or vector of retention time indices
#      rt.observed: A matrix or vector of observed retention times.
#      ri.standard: A vector of retention index standard.
# notes:
#   1. x and rt.observed must have the same number of columns. It is assumed that
#      there is a correspondence between columns.
#   2. the number of rows of rt.observed must be equal to length of ri.standard.

ri2rt <- function(riTime, rt.observed, ri.standard) {
	if(is.null(dim(riTime))) {
		stopifnot(length(rt.observed) == length(ri.standard))
		return(rt2ri(riTime, ri.standard, rt.observed))
	}
	stopifnot(!is.null(dim(rt.observed)))
	stopifnot(ncol(riTime) == ncol(rt.observed))
  stopifnot(length(ri.standard) == nrow(rt.observed))
  xo <- sapply(1:ncol(riTime), function(n) rt2ri(riTime[,n], ri.standard, rt.observed[,n]))
  if(!is.null(colnames(riTime))) colnames(xo) <- colnames(riTime)
  if(!is.null(rownames(riTime))) rownames(xo) <- rownames(riTime)
	return(xo)
}
