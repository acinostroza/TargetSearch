# tests for class tsSample

# make an object
cdf <- sprintf("file_%d.nc4", 1:6)
data <- data.frame(CDF=cdf, Treatment=c('A','A','A', 'B','B','B'))
smp <- new('tsSample', CDFfiles=cdf, RIpath='.', days=1, data=data)
expect_true(validObject(smp))
expect_equal(length(smp), 6L)
expect_equal(CDFfiles(smp), cdf)

# make another
cdf <- sprintf("file_%d.nc4", 7:12)
data <- data.frame(CDF=cdf, Treatment=c('C','C','C', 'D','D','D'))
smp2 <- new('tsSample', CDFfiles=cdf, RIpath='.', days=2, data=data)
expect_true(validObject(smp2))
expect_equal(length(smp2), 6L)
expect_equal(CDFfiles(smp2), cdf)

# test combination
s <- c(smp, smp2)
expect_true(validObject(s))
expect_equal(length(s), 12L)
expect_equal(CDFfiles(s), sprintf("file_%d.nc4", seq(12)))
