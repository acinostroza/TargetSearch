# test text2bin and bin2text functions

require(TargetSearchData)
tmp <- tempdir()

in.files <- dir(file.path(find.package("TargetSearchData"), "gc-ms-data"), pattern="RI_", full=TRUE)
in.files <- sample(in.files, 3)
expect_equal(file.exists(in.files), c(TRUE, TRUE, TRUE))

out.files <- sub(".txt", ".dat", basename(in.files))
test <- text2bin(in.files, file.path(tmp, out.files))
expect_equal(out.files, basename(test))

test2 <- bin2text(test)
expect_equal(basename(in.files), basename(test2))

# converted files should be equal to original
a <- lapply(in.files, read.delim, nrows=10)
b <- lapply(test2, read.delim, nrows=10)

expect_equal(a, b)
