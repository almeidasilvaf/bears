
sinfo <- structure(
    list(BioSample = c("SAMN03000709", "SAMN02047242", 
                       "SAMN02141349", "SAMN01779511", "SAMN01779511", "SAMN01779510", 
                       "SAMN01779510", "SAMN01779510", "SAMN01779509", "SAMN03003058", 
                       "SAMN03169152"), 
         Experiment = c("SRX684378", "SRX268041", "SRX276001", 
                        "SRX200138", "SRX200138", "SRX200137", "SRX200137", "SRX200137", 
                        "SRX200136", "SRX1880062JGI-SRA-14008", "SRX763923"), 
         Run = c("SRR1555233", 
                 "SRR830199", "SRR847642", "SRR605677", "SRR605678", "SRR605674", 
                 "SRR605675", "SRR605676", "SRR605671", "SRR3722300", "SRR1657460"
         ), 
         BioProject = c("PRJNA258604", "PRJNA197379", "PRJNA178155", 
                        "PRJNA178155", "PRJNA178155", "PRJNA178155", "PRJNA178155", "PRJNA178155", 
                        "PRJNA178155", "PRJNA250761", "PRJNA266538"), 
         Instrument = c("Illumina HiSeq 2000", 
                        "Illumina HiSeq 2000", "Illumina HiSeq 2000", "Illumina HiSeq 2000", 
                        "Illumina HiSeq 2000", "Illumina HiSeq 2000", "Illumina HiSeq 2000", 
                        "Illumina HiSeq 2000", "Illumina HiSeq 2000", "Illumina HiSeq 2000", 
                        "Illumina HiSeq 2000"), 
         Layout = c("SINGLE", "PAIRED", "SINGLE", 
                    "SINGLE", "SINGLE", "SINGLE", "SINGLE", "SINGLE", "SINGLE", "PAIRED", 
                    "SINGLE"), 
         Selection_method = c("cDNA", "cDNA", "RANDOM", "RANDOM", 
                              "RANDOM", "RANDOM", "RANDOM", "RANDOM", "RANDOM", "RT-PCR", "cDNA"
         )), 
    row.names = c(462L, 763L, 1202L, 
                  1255L, 1256L, 1258L, 1259L, 1260L, 1262L, 1298L, 1318L), 
    class = "data.frame"
)

test_that("create_sample_info() creates a data frame of search results", {
    term <- "PRJNA80173[GPRJ]"
    df <- create_sample_info(term)
    expect_equal(nrow(df), 4)
    expect_equal(ncol(df), 17)
    expect_equal(class(df), "data.frame")
})

test_that("get_url_ena() correcly gets ENA's URL to FASTQ files", {
    idx <- sample(seq_len(nrow(sinfo)), 1)
    urls <- get_url_ena(sinfo[idx, ])
    expect_equal(class(urls), "list")
})
