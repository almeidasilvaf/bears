
#----Load data------------------------------------------------------------------
data(sample_info)

## Update PATH
Sys.setenv(
    PATH = paste(
        Sys.getenv("PATH"),
        "/opt/STAR-2.7.9a/bin/Linux_x86_64_static/",
        "/opt/bin/",
        "/opt/salmon-1.5.2_linux_x86_64/bin/",
        "/opt/kallisto/",
        "/opt/subread-2.0.3-Linux-x86_64/bin/",
        "/opt/stringtie-2.1.7.Linux_x86_64/",
        "/opt/taco-v0.7.3.Linux_x86_64/",
        sep = ":"
    )
)


#----Start tests----------------------------------------------------------------
test_that("get_fastq_paths() returns paths to FASTQ files", {
    
    fastqdir <- system.file("extdata", package = "bears")
    
    g1 <- get_fastq_paths(fastqdir, "SRR1039508", cmd = "P1")
    g2 <- get_fastq_paths(fastqdir, "SRR1039508", cmd = "P2")
    g3 <- get_fastq_paths(fastqdir, "SRR1039508", cmd = "S")
    g4 <- get_fastq_paths(fastqdir, "SRR1039509", cmd = "P1")
    g5 <- get_fastq_paths(fastqdir, "SRR1039509", cmd = "P2")
    
    expect_equal(class(g1), "character")    
    expect_equal(class(g2), "character")
    expect_null(g3)
    expect_null(g4)
    expect_null(g5)
})



