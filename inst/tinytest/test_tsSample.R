# tests for class tsSample

# make an object
cdf <- sprintf("file_%d.nc4", 1:6)
data <- data.frame(CDF=cdf, Treatment=c('A','A','A', 'B','B','B'))
expect_silent(smp <- new('tsSample', CDFfiles=cdf, RIpath='.', days=1, data=data))
expect_true(validObject(smp))
expect_equal(length(smp), 6L)
expect_equal(CDFfiles(smp), cdf)

# make another
cdf <- sprintf("file_%d.nc4", 7:12)
data <- data.frame(CDF=cdf, Treatment=c('C','C','C', 'D','D','D'))
expect_silent(smp2 <- new('tsSample', CDFfiles=cdf, RIpath='.', days=2, data=data))
expect_true(validObject(smp2))
expect_equal(length(smp2), 6L)
expect_equal(CDFfiles(smp2), cdf)

# test combination
expect_silent(s <- c(smp, smp2))
expect_true(validObject(s))
expect_equal(length(s), 12L)
expect_equal(CDFfiles(s), sprintf("file_%d.nc4", seq(12)))

# test combination unequal columns
cdf <- sprintf("file_%d.nc4", 13:16)
data <- data.frame(CDF=cdf, Treatment=c('E','E', 'F','F'), TimePoint=1:4)
expect_silent(smp3 <- new('tsSample', CDFfiles=cdf, RIpath='.', days=2, data=data))
expect_true(validObject(smp3))

expect_silent(s <- c(smp, smp2, smp3))
expect_true(validObject(s))
expect_equal(length(s), 16L)
expect_equal(CDFfiles(s), sprintf("file_%d.nc4", seq(16)))
expect_equal(ncol(sampleData(s)), 3L)
