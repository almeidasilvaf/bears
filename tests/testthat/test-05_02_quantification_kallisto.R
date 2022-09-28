
#----Load data------------------------------------------------------------------
data(tx2gene)
data(sample_info)
sample_info <- rbind(
    sample_info, sample_info
)
sample_info$BioSample[2] <- "SAMN02422670"
sample_info$Experiment[2] <- "SRX384346"
sample_info$BioProject[2] <- "PRJNA229999"
sample_info$Run[2] <- "SRR1039509"
sample_info$Layout[2] <- "SINGLE" 

qc_table <- summary_stats_fastp(system.file("extdata", package = "bears"))
qc_table$Sample <- "SRR1039509"

## Create fake single-end FASTQ file
fake_fastq <- c(
    "@SRR1039508.2486",
    "TTCAACCTTTCATGTAACAAAACTTATAAAATTCCTTTAGCTCTTCCTTTTTCTAAAGAAAAA",
    "+",
    "HJJJJJJJJJJJJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJJJJJJ",
    "@SRR1039508.14392",
    "CTGGCACACCTGATCAAATTTCTCCTCCATCAAGTCTCTGCAACACCGGCTCTCCACCAGGGT",
    "+",
    "HJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIJGIJJHHHHFFFFFDDD;",
    "@SRR1039508.14861",
    "GACATTACTTTCACTTGTGTTTCATATTCCGTAAAACGAACAAAGCCAAACCCCTTTGAATGA",
    "+",
    "HJJJJJJJJJJJJJJJJHHIJJJJJJJJJIJJIJIJJJJJJJJJJJJJJJIAHHHHFFFFFFF",
    "@SRR1039508.18596",
    "AGCGGCCAGATCGGCTTCATGCTGTGCAGCAAGAACCCGAGCACGAACTTCCAGGAGCCGGTG",
    "+",
    "HJJJJJJJJJJJJJJJJI>FGGGEBFB<CF@HGFHADHEEEFDED@BB<C@CCCDBBBD9@;@",
    "@SRR1039508.25193",
    "TCATACACTGAGACACTGGATTCCTGTAGCGAAGCCCACACGCCCCTGGAAACTGGGCTGTAA",
    "+"
)

fastqdir <- file.path(tempdir(), "fastqdir")
if(!dir.exists(fastqdir)) { dir.create(fastqdir, recursive = TRUE) }

## Write fake FASTQ to `fastqdir`
gz1 <- gzfile(file.path(fastqdir, "SRR1039509.fastq.gz"), "w")
writeLines(text = fake_fastq, con = gz1)
close(gz1)

## Copy FASTQ files in inst/extdata to `fastqdir`
fqfiles <- list.files(
    system.file("extdata", package = "bears"), pattern = ".fastq.gz",
    full.names = TRUE
)
c1 <- file.copy(from = fqfiles[1], to = fastqdir)
c2 <- file.copy(from = fqfiles[2], to = fastqdir)

# Path to files and directories used in this set of tests
kallistoindex <- file.path(tempdir(), "kallistoidx")
transcriptome_path <- system.file(
    "extdata", "Homo_sapiens.GRCh37.75_subset_transcripts.fa.gz", 
    package = "bears"
)
kallistodir <- file.path(tempdir(), "kallisto_quant")

## Create a random SummarizedExperiment object
nrows <- 200
ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colData <- data.frame(
    Treatment = rep(c("ChIP", "Input"), 3),
    row.names = LETTERS[1:6]
)
se_random <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), colData = colData
)

## Change PATH
Sys.setenv(
    PATH = paste(
        Sys.getenv("PATH"),
        paste0(Sys.getenv("HOME"), "/.local/bin"),
        paste0(
            Sys.getenv("HOME"), 
            "/Documents/Programs/salmon-1.9.0_linux_x86_64/bin/"
        ),
        "/opt/STAR-2.7.9a/bin/Linux_x86_64_static/",
        "/opt/salmon-1.9.0_linux_x86_64/bin/",
        "/opt/bin/",
        "/opt/kallisto/",
        "/opt/subread-2.0.3-Linux-x86_64/bin/",
        "/opt/stringtie-2.1.7.Linux_x86_64/",
        "/opt/taco-v0.7.3.Linux_x86_64/",
        sep = ":"
    )
)


#----Start tests----------------------------------------------------------------
test_that("kallisto_index(), kallisto_quantify(), and kallisto2se() work", {
    
    ki <- data.frame()
    kq <- data.frame()
    se_gene <- se_random
    se_tx <- se_random
    se_both <- list()
    
    if(salmon_is_installed()) {
        ki <- kallisto_index(kallistoindex, transcriptome_path)
        kq <- kallisto_quantify(
            sample_info, 
            qc_table = qc_table, 
            filtdir = fastqdir, kallistoindex = kallistoindex, 
            kallistodir = kallistodir
        )
        d <- fs::dir_delete(file.path(kallistodir, "SAMN02422670"))
        se_gene <- kallisto2se(sample_info[1, ], "gene", kallistodir, tx2gene)
        se_tx <- kallisto2se(sample_info[1, ], "transcript", kallistodir, tx2gene)
        se_both <- kallisto2se(sample_info[1, ], "both", kallistodir, tx2gene)
        
        expect_error(
            kallisto2se(sample_info[1, ], "error", kallistodir, tx2gene)
        )
    }
    
    expect_equal(class(ki), "data.frame")
    expect_equal(class(kq), "data.frame")
    
    expect_true("SummarizedExperiment" %in% class(se_gene))
    expect_true("SummarizedExperiment" %in% class(se_tx))
    expect_equal(class(se_both), "list")
})


test_that("run2biosample_kallisto() handles technical replicates", {
    
    r1 <- run2biosample_salmon(sample_info, fastqdir)
    
    expect_equal(class(r1), "list")
    expect_equal(length(r1), 2)
    expect_equal(names(r1), c("paired", "single"))
})


