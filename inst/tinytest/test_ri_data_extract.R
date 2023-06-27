require('TargetSearchData')

files <- tsd_rifiles()

expect_true(all(file.exists(files)))

# pure R implementation that works exactly as ri_data_extract
ri_data_extract_text <- TargetSearch:::ri_data_extract_text

a <- lapply(files, ri_data_extract,  204, c(200, 300), useRT=TRUE)
b <- lapply(files, ri_data_extract_text,  204, c(200, 300), useRT=TRUE)
expect_identical(a, b)

a <- lapply(files, ri_data_extract, 116, c(200, 220) * 1000)
b <- lapply(files, ri_data_extract_text, 116, c(200, 220) * 1000)
expect_identical(a, b)

a <- lapply(files, ri_data_extract, 1234, c(200, 220) * 1000)
expect_equal(sapply(a, nrow), integer(length(files)))

# test mock data
source('mock_functions.R')
data <- mock_peak_data(mz_range=c(100,150), time_range=c(300, 400), n_peaks=200, n_scans=250, ri_factor=1004)
myfile <- mock_write_file(data)
txtfile <- bin2text(myfile, paste0(myfile, ".txt"))

a <- ri_data_extract(myfile, c(100:150), c(299, 401), useRT=TRUE)
b <- ri_data_extract_text(txtfile, c(100:150), c(299, 401), useRT=TRUE)
expect_equivalent(a, b)

a <- ri_data_extract(myfile, c(100:150), c(299, 401) * 1004)
b <- ri_data_extract_text(txtfile, c(100:150), c(299, 401) * 1004)
expect_equivalent(a, b)
