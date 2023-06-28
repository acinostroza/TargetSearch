if ( requireNamespace("tinytest", quietly=TRUE) && requireNamespace("TargetSearchData", quietly=TRUE)) {
    tinytest::test_package("TargetSearch", color=FALSE, verbose=1)
}
