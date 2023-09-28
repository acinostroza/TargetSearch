# test internal function writeRIFile

source('mock_functions.R')

writeRIFile <- TargetSearch:::writeRIFile

data <- mock_peak_data(c(100,200), c(300, 400), 500, 401)
ritxt <- tempfile(fileext='.txt')
ribin <- tempfile(fileext='.bin')

expect_true(writeRIFile(ritxt, data, data$Index, data$massRange, 'text'))
expect_true(writeRIFile(ribin, data, data$Index, data$massRange, 'bin'))
expect_error(writeRIFile(tempfile(), data, data$Index, data$massRange, 'txt'))

mass <- sample(seq(data$massRange[1], data$massRange[2]), 20)
time <- sample(data$Time, 1) + c(-20, 20)
a <- ri_data_extract(ritxt, mass, time, useRT=TRUE)
b <- ri_data_extract(ribin, mass, time, useRT=TRUE)
expect_equivalent(a, b)

# check file signatures
expect_identical(readLines(ritxt, 1L), TargetSearch:::get.file.header())
sig <- as.raw(c(42, 243, 27, 10, 183, 75, 1, 0))
expect_identical(readBin(ribin, 'raw', 8L), sig)
