
#----Load data------------------------------------------------------------------
data(sample_info)

## Update PATH
Sys.setenv(
    PATH = paste(
        Sys.getenv("PATH"),
        "/opt/FastQC/fastqc",
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

test_that("get_fastqc_paths() returns paths to FASTQC files", {
    
    fastqdir <- system.file("extdata", package = "bears")
    fastqcdir <- system.file("extdata", package = "bears")
    
    g1 <- get_fastqc_paths(fastqdir, fastqcdir, "SRR1039508", cmd = "P1")
    g2 <- get_fastqc_paths(fastqdir, fastqcdir, "SRR1039508", cmd = "P2")
    g3 <- get_fastqc_paths(fastqdir, fastqcdir, "SRR1039508", cmd = "S")
    g4 <- get_fastqc_paths(fastqdir, fastqcdir, "SRR1039509", cmd = "P1")
    g5 <- get_fastqc_paths(fastqdir, fastqcdir, "SRR1039509", cmd = "P2")
    
    expect_equal(class(g1), "list")    
    expect_equal(class(g2), "list")
    expect_equal(class(g3), "list")
    expect_equal(class(g4), "list")
    expect_equal(class(g5), "list")
})

test_that("run_fastqc() returns a 2-column status data frame", {
    
    fastqdir <- system.file("extdata", package = "bears")
    fastqcdir <- tempdir()
    
    f <- data.frame()
    if(fastqc_is_installed()) {
        f <- run_fastqc(sample_info[1, ], fastqdir, fastqcdir)
    }

    expect_equal(class(f), "data.frame")
})


test_that("multiqc() returns a QC summary data frame", {
    
    m1 <- data.frame()
    m2 <- data.frame()
    
    idir <- system.file("extdata", package = "bears")
    
    if(multiqc_is_installed()) {
        m1 <- multiqc(
            idir,
            outdir = file.path(tempdir(), "fastq_cout"),
            runon = "fastqc"
        )
        
        m2 <- multiqc(
            idir,
            outdir = file.path(tempdir(), "star_cout"),
            runon = "star"
        )
    }
    
    expect_equal(class(m1), "data.frame")
    expect_equal(class(m2), "data.frame")
})


