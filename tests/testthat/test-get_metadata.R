
#---Load example data-----------------------------------------------------------
data(sample_info)

#----Start tests----------------------------------------------------------------
test_that("create_sample_info() creates a data frame of search results", {
    term <- "PRJNA80173[GPRJ]"
    df <- create_sample_info(term)
    expect_equal(nrow(df), 4)
    expect_equal(ncol(df), 17)
    expect_equal(class(df), "data.frame")
})

test_that("get_read_count() returns a data frame of expected read number", {
    readcount <- get_read_count(sample_info, sample_info$Run)
    
    expect_equal(class(readcount), "data.frame")
    expect_equal(ncol(readcount), 2)
    expect_equal(class(readcount$Run), "character")
    expect_equal(class(readcount$Reads), "numeric")
})
