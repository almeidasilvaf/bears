
#----Create example data--------------------------------------------------------
sample_metadata <- data.frame(
    BioSample = c(
        "SAMN03000709", "SAMN02047242", "SAMN02141349", "SAMN01779511"
    ), 
    Experiment = c("SRX684378", "SRX268041", "SRX276001", "SRX200138"),
    Run = c("SRR1555233", "SRR830199", "SRR847642", "SRR605677"), 
    BioProject = c("PRJNA258604", "PRJNA197379", "PRJNA178155", "PRJNA178155"),
    Instrument = c(rep("Illumina HiSeq 2000", 4)), 
    Layout = c("SINGLE", "PAIRED", "SINGLE", "SINGLE"), 
    Selection_method = c("cDNA", "cDNA", "RANDOM", "RANDOM")
)

fastqdir <- system.file("extdata", package = "bears") 
read_count <- data.frame(Run = "SRR1039508", Reads = 7097)


#----Start tests----------------------------------------------------------------
test_that("get_url_ena() correcly gets ENA's URL to FASTQ files", {
    idx <- sample(seq_len(nrow(sample_metadata)), 1)
    urls <- get_url_ena(sample_metadata[idx, ])
    expect_equal(class(urls), "character")
})

test_that("check_downloads() returns a data frame with download status", {
    dcheck <- check_downloads(sample_info, fastqdir, read_count)
    
    expect_equal(class(dcheck), "data.frame")
    expect_equal(ncol(dcheck), 3)
    expect_equal(names(dcheck), c("Run", "Status", "Issue"))
})
