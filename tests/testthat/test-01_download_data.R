
options(timeout = 60000)

#----Create example data--------------------------------------------------------
sample_metadata <- data.frame(
    BioSample = c(
        "SAMN03000709", "SAMN02047242", "SAMN02141349", "SAMN01779511",
        "SAMN000001"
    ), 
    Experiment = c(
        "SRX684378", "SRX268041", "SRX276001", "SRX200138", "SRX0001"
    ),
    Run = c(
        "SRR1555233", "SRR830199", "SRR847642", "SRR605677", "SRR1039508"
    ), 
    BioProject = c(
        "PRJNA258604", "PRJNA197379", "PRJNA178155", "PRJNA178155", 
        "PRJNA178155"
    ),
    Instrument = c(rep("Illumina HiSeq 2000", 5)), 
    Layout = c("SINGLE", "PAIRED", "SINGLE", "SINGLE", "PAIRED"), 
    Selection_method = c("cDNA", "cDNA", "RANDOM", "RANDOM", "cDNA")
)

fastqdir <- system.file("extdata", package = "bears") 


#----Start tests----------------------------------------------------------------
test_that("fastq_exists() checks if FASTQ files exist in fastqdir", {

    f1 <- fastq_exists(sample_metadata, fastqdir)
    f2 <- fastq_exists(sample_metadata, fastqdir, collapse_pe = FALSE)
    
    expect_equal(class(f1), "data.frame")
    expect_equal(class(f2), "data.frame")
})


test_that("get_url_ena_api() works", {
    
    fake_sample_info <- data.frame(Run = "fakerun")
    g <- get_url_ena_api(fake_sample_info)
    expect_null(g)
})


test_that("get_url_ena_iterative() works", {
    meta2 <- sample_metadata
    meta2$Run[2] <- "ERR550666"
    meta2$Run[3] <- "SRR000000"
    
    f1 <- get_url_ena_iterative(sample_metadata[c(1,5), ])
    f2 <- get_url_ena_iterative(meta2[2, ])
    f3 <- get_url_ena(meta2[3, ])
    
    expect_equal(class(f1), "character")
    expect_equal(class(f2), "character")
    expect_null(f3)
})


test_that("get_url_ena() correcly gets ENA's URL to FASTQ files", {
    idx <- sample(seq_len(nrow(sample_metadata)), 1)

    urls <- get_url_ena(sample_metadata[idx, ])
    urls2 <- get_url_ena(sample_metadata[idx, ], link_from = "iterative")
    
    expect_error(get_url_ena(sample_metadata[idx, ], link_from = "error"))
    expect_equal(class(urls), "character")
    expect_equal(class(urls2), "character")
})


test_that("download_from_ena() downloads FASTQ files", {
    
    sample_info <- data.frame(
        BioSample = "SAMN01924555",
        Experiment = "SRX245306",
        Run = "SRR926397",
        BioProject = "PRJNA190191", Instrument = "Illumina HiSeq 2000", 
        Layout = "PAIRED"
    )
    fastqdir <- tempdir()
    
    d <- download_from_ena(
        sample_info, fastqdir = fastqdir
    )
    
    expect_equal(class(d), "data.frame")
})

test_that("check_md5() checks file integrity and returns a data frame", {
    
    c1 <- check_md5("SRR926397", tempdir())
    
    expect_equal(class(c1), "data.frame")
})


