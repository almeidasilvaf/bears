test_that("create_sample_info() creates a data frame of search results", {
    term <- "PRJNA80173[GPRJ]"
    df <- create_sample_info(term)
    expect_equal(nrow(df), 4)
    expect_equal(ncol(df), 17)
    expect_equal(class(df), "data.frame")
})
