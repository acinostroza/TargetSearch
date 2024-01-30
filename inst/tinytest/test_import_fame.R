# test ImportFameSettings

# random RI standard definition
mass <- 85
std  <- c(123, 456, 789)
lim  <- cbind(a=c(10, 20, 30), b=c(20, 30, 40))

rimfile <- tempfile()
write.table(cbind(lim, std), file=rimfile, sep="\t")

expect_silent(x <- ImportFameSettings(rimfile))
expect_silent(y <- ImportFameSettings(rimfile, mass=mass))

expect_equal(rimMass(x), 87)
expect_equal(rimMass(y), mass)
expect_equal(rimStandard(x), rimStandard(y))
expect_equal(rimLimits(x), rimLimits(y))

write.table(cbind(lim, std, mass), file=rimfile, sep="\t")

expect_silent(x <- ImportFameSettings(rimfile))
expect_equal(rimMass(x),rep(mass, length(std)))
