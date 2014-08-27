## unit tests will not be done if testthat is not available
if(require("testthat", quietly=TRUE)) {
    library(hiReadsProcessor)
    data(FLX_seqProps)
    test_package("hiReadsProcessor")    
} else {
    warning("cannot run unit tests -- package 'testthat' is not available")
}
