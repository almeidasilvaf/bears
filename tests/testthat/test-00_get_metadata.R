
#---Load example data-----------------------------------------------------------
data(sample_info)

#----Start tests----------------------------------------------------------------
test_that("check_empty() works", {
    c1 <- check_empty(c("test1", "test2"))
    c2 <- check_empty(character(0))
    c3 <- check_empty(list("test1", "test2"))
    
    expect_equal(c1, "test1")
    expect_true(is.na(c2))
    expect_true(is.na(c3))
})

test_that("create_meta_df() create a data frame of metadata", {
    
    fake_list <- list(
        biosample = "bla", experiment = "bla", run = "bla", 
        tissue = "bla", pubmed = "bla", bioproject = "bla", 
        instrument = "bla", layout = "bla", selection = "bla", 
        srasample = "bla", srastudy = "bla", treatment = "bla",
        cultivar = "bla", title = "bla", abstract = "bla", date = "bla",
        origin = "bla"
    )
    fake_list2 <- fake_list
    fake_list2$run <- NULL
    
    c1 <- create_meta_df(fake_list)
    c2 <- create_meta_df(fake_list2)
    
    expect_equal(class(c1), "data.frame")
    expect_true(is.na(c2$Run))
})

test_that("sra_xml2df() creates a data frame from XML output", {
    
    s1 <- sra_xml2df(id = "551386")
    expect_equal(class(s1), "data.frame")
})


test_that("create_sample_info() creates a data frame of search results", {
    term <- "PRJNA80173[GPRJ]"
    df <- create_sample_info(term)
    expect_equal(nrow(df), 4)
    expect_equal(ncol(df), 17)
    expect_equal(class(df), "data.frame")
})
