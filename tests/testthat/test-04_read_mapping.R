

#----Load data------------------------------------------------------------------
data(sample_info)
sample_info <- rbind(
    sample_info, sample_info
)
sample_info$BioSample[2] <- "SAMN02422670"
sample_info$Experiment[2] <- "SRX384346"
sample_info$BioProject[2] <- "PRJNA229999"
sample_info$Run[2] <- "SRR1039509"
sample_info$Layout[2] <- "SINGLE" 

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
        "/opt/bin/",
        "/opt/kallisto/",
        "/opt/subread-2.0.3-Linux-x86_64/bin/",
        "/opt/stringtie-2.1.7.Linux_x86_64/",
        "/opt/taco-v0.7.3.Linux_x86_64/",
        sep = ":"
    )
)

## Path to files and directories used in this set of tests

genome_path <- system.file(
    "extdata", "Homo_sapiens.GRCh37.75_subset.fa", package = "bears"
)
gff_path <- system.file(
    "extdata", "Homo_sapiens.GRCh37.75_subset.gtf", package = "bears"
)
mappingdir <- file.path(tempdir(), "mappingdir")
qc_table <- summary_stats_fastp(system.file("extdata", package = "bears"))


#----Start tests----------------------------------------------------------------
test_that("star_genome_index() and star_align() work", {
    
    g <- data.frame()
    a <- data.frame()
    
    if(star_is_installed()) {
        g <- star_genome_index(genome_path, gff_path, mappingdir)
        a <- star_align(
            sample_info[1,], fastqdir, qc_table, mappingdir, gff_path
        )
    }
    
    expect_equal(class(g), "data.frame")
    expect_equal(class(a), "data.frame")
})

test_that("star_reads() creates a list of paths for STAR", {
    
    s <- star_reads(sample_info, fastqdir)
    
    expect_equal(class(s), "list")
    expect_equal(length(s), 2)
})

test_that("mapping_pass() filters only samples that pass minimum criteria", {
    
    mapping_qc <- summary_stats_star(system.file("extdata", package = "bears"))
    m <- mapping_pass(mapping_qc, sample_info[1, ])
    
    expect_equal(class(m), "data.frame")
    expect_equal(nrow(m), 1)
})

