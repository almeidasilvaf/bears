data(sample_info)

test_that("run_fastqc works", {
    fq <- system.file("extdata", package="bear")
    fqc <- tempdir()
    n <- run_fastqc(sample_info[1, ], fastqdir = fq, fastqcdir = fqc)
    expect_equal(n, NULL)
})
