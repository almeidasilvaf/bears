
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

## Create dbdir for SortMeRNA
dbdir <- file.path(tempdir(), "rrna")
dir.create(dbdir)
rrna_file <- system.file("extdata", "bac_16s_subset.fa", package = "bears")
copy <- file.copy(from = rrna_file, to = dbdir)

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
test_that("trim_reads() works and returns a 2-column status data frame", {
    
    filtdir <- file.path(tempdir(), "filtdir")
    qcdir <- file.path(tempdir(), "qcdir")
    
    t1 <- data.frame()

    if(fastp_is_installed()) {
        t1 <- trim_reads(
            sample_info, fastqdir, filtdir, qcdir, delete_raw = FALSE
        )
    }
    
    expect_equal(class(t1), "data.frame")
})


test_that("rrna_ref_args() creates a vector of SortMeRNA args", {
    
    dbdir <- system.file("extdata", package = "bears")
    r1 <- rrna_ref_args(dbdir)
    
    expect_equal(class(r1), "character")
})

test_that("clean_sortmerna() cleans result directory", {
    
    filtdir <- file.path(tempdir(), "sim_rrnadir")
    if(!dir.exists(filtdir)) { dir.create(filtdir, recursive = TRUE) }
    
    c <- clean_sortmerna(filtdir)
    expect_null(c)
})

test_that("remove_rrna() removes rRNA from FASTQ files", {
    
    r1 <- data.frame()
    
    if(sortmerna_is_installed()) {
        r1 <- remove_rrna(
            sample_info, 
            fastqdir = fastqdir,
            filtdir = file.path(tempdir(), "filtdir"),
            rrna_db_dir = dbdir
        )
    }
    
    expect_equal(class(r1), "data.frame")
})

