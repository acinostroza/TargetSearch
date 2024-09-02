# tests for class tsRim

# test creation and length method
expect_silent(rim <- new('tsRim', limits=cbind(1:6, 11:16), standard=100 * 1:6, mass=87))
expect_true(validObject(rim))
expect_equal(length(rim), 6L)

# test combination with `c`
expect_silent(rim2 <- new('tsRim', limits=cbind(rnorm(7), rnorm(7) + 5), standard=100 * runif(7), mass=3))
expect_true(validObject(rim2))
expect_equal(length(rim2), 7L)

expect_silent(x <- c(rim, rim2))
expect_equal(length(x), 13L)
expect_true(validObject(x))
