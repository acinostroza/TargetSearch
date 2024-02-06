#
# check importing different type of files
#

# functions:
## write file
wrf <- function(x, f) write.table(x, file=f, sep="\t", quote=FALSE)
## format spectrum data
frm <- Vectorize(function(x) {
    if(is.matrix(x)) paste(sprintf("%d:%d", x[,1], x[,2]), collapse=" ") else "NA"
}, USE.NAMES=FALSE)

libfile <- tempfile()

# test1: no mass (should fail)
libdata <- data.frame(Name=sprintf("met_%d", 1:2), RI=c(100, 200))
wrf(libdata, libfile)
expect_error(lib <- ImportLibrary(libfile))

# test2: minimal library
libdata <- data.frame(Name=sprintf("met_%d", 1:2), RI=c(100, 200), SEL_MASS=c(117, 147))
wrf(libdata, libfile)
expect_silent(lib <- ImportLibrary(libfile))
expect_equivalent(libName(lib), libdata$Name)
expect_equivalent(libRI(lib), libdata$RI)

# test3: minimal library with top mass
libdata <- data.frame(Name=sprintf("met_%d", 1:2), RI=c(100, 200), TOP_MASS=c(117, 147))
wrf(libdata, libfile)
expect_silent(lib <- ImportLibrary(libfile))

# test4: minimal example with spectrum
sp <- list(cbind(c(147, 167, 201), c(999, 10, 235)), cbind(c(234, 300), c(999, 10)))
libdata <- data.frame(Name=sprintf("met_%d", 1:2), RI=c(100, 200), SPECTRUM=frm(sp))
wrf(libdata, libfile)
expect_silent(lib <- ImportLibrary(libfile))
expect_equivalent(libName(lib), libdata$Name)
expect_equivalent(libRI(lib), libdata$RI)
expect_equivalent(spectra(lib), sp)

# test5: invalid spectrum data (should fail)
sp <- list(cbind(c(147, 167, 201), c(999, 10, 235)), numeric(0))
libdata <- data.frame(Name=sprintf("met_%d", 1:2), RI=c(100, 200), SPECTRUM=frm(sp))
wrf(libdata, libfile)
expect_error(lib <- ImportLibrary(libfile))

# test6: add selective masses so it passes
libdata$SEL_MASS <- c(117, 147)
wrf(libdata, libfile)
expect_silent(lib <- ImportLibrary(libfile))
expect_equivalent(spectra(lib), sp)

# test7: add quan mass
libdata$QUANT_MASS <- c(143, 201)
wrf(libdata, libfile)
expect_silent(lib <- ImportLibrary(libfile))
expect_equivalent(quantMass(lib), libdata$QUANT_MASS)


# check the objects in example data can be generated identically
require(TargetSearchData)
data(TSExample)
lib <- ImportLibrary(tsd_file_path('library.txt'))
expect_identical(lib, refLibrary)
