
#----Load data------------------------------------------------------------------


#----Start tests----------------------------------------------------------------
test_that("is_valid() works", {
    v <- is_valid(cmd = "echo", args = c("Hello"))
    
    expect_equal(class(v), "logical")
})


test_that("*_is_installed() functions return a logical scalar", {

    c1 <- fastp_is_installed()    
    c2 <- multiqc_is_installed()
    c4 <- sortmerna_is_installed()
    c5 <- star_is_installed()
    c6 <- salmon_is_installed()
    c7 <- kallisto_is_installed()
    c8 <- subread_is_installed()
    c9 <- taco_is_installed()
    c10 <- rseqc_is_installed()
    c11 <- stringtie_is_installed()
    
    expect_equal(class(c1), "logical")
    expect_equal(class(c2), "logical")
    expect_equal(class(c4), "logical")
    expect_equal(class(c5), "logical")
    expect_equal(class(c6), "logical")
    expect_equal(class(c7), "logical")
    expect_equal(class(c8), "logical")
    expect_equal(class(c9), "logical")
    expect_equal(class(c10), "logical")
    expect_equal(class(c11), "logical")
})
