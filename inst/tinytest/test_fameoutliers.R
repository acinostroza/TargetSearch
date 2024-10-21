# test Fame outliers

# make sample object
cdf <- sprintf("file_%d.nc4", 1:25)
data <- data.frame(CDF=cdf, Treatment=rep(LETTERS[1:5], each=5))
days <- data$Treatment
expect_silent(smp <- new('tsSample', CDFfiles=cdf, RIpath='.', days=days, data=data))

# make RImatrix
RT <- c(252, 311, 369, 432, 554, 638, 715) # sample RT
RImatrix <- outer(RT, rnorm(length(smp)), '+')
pdffile <- tempfile(fileext='.pdf')

# simple tests in which we test start/end days and default behavior
expect_silent(out <- FAMEoutliers(smp, RImatrix, pdffile))
expect_silent(out <- FAMEoutliers(smp, RImatrix, pdffile, startDay=NA, endDay=NA))
expect_silent(out <- FAMEoutliers(smp, RImatrix, pdffile, startDay=c('A'), endDay=c('E')))
expect_silent(out <- FAMEoutliers(smp, RImatrix, pdffile, startDay=c('A', 'C'), endDay=c('D', 'E')))
