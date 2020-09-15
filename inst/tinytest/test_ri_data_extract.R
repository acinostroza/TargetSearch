require('TargetSearchData')

path <- file.path(find.package("TargetSearchData"), "gc-ms-data")
id <- c(4,6,7,8,9,11,12,15)

files <- file.path(path, sprintf('RI_7235eg%02d.txt', id))

expect_equal(file.exists(files), rep(TRUE, 8))

a <- lapply(files, ri_data_extract,  204, c(200, 300), useRT=TRUE)
b <- c(26L, 18L, 22L, 28L, 23L, 33L, 27L, 20L)
expect_equal(sapply(a, nrow), b)

a <- lapply(files, ri_data_extract, 116, c(200, 220) * 1000)
b <- c(1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L)
expect_equal(sapply(a, nrow), b)

a <- lapply(files, ri_data_extract, 1234, c(200, 220) * 1000)
expect_equal(sapply(a, nrow), integer(8))

# test mock data
source('mock_functions.R')
data <- mock_peak_data(mz_range=c(100,150), time_range=c(300, 400), n_peaks=200, n_scans=250, ri_factor=1004)
myfile <- mock_write_file(data)
m <- which(data$Peaks > 0)
ret <- ri_data_extract(myfile, c(100:150), c(299, 401), useRT=TRUE)
expect_equal( ret[, 'Intensity'] ,  data$Peaks[ m ] )
ret <- ri_data_extract(myfile, c(100:150), c(299, 401) * 1004)
expect_equal( ret[, 'Intensity'] ,  data$Peaks[ m ] )
