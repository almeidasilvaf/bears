

#----Load data------------------------------------------------------------------
data(sample_info)
data(fastqc_table)

## Create dbdir for SortMeRNA
dbdir <- file.path(tempdir(), "rrna")
dir.create(dbdir)
rrna_file <- system.file("extdata", "bac_16s_subset.fa", package = "bears")
copy <- file.copy(from = rrna_file, to = dbdir)

#----Start tests----------------------------------------------------------------
test_that("trim_reads() works and returns a 2-column status data frame", {
    
    fastqdir <- system.file("extdata", package = "bears")
    filtdir <- file.path(tempdir(), "filtdir")
    
    t1 <- data.frame()

    if(trimmomatic_is_installed()) {
        t1 <- trim_reads(
            sample_info, fastqc_table, fastqdir, filtdir
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
            fastqdir = system.file("extdata", package = "bears"),
            filtdir = file.path(tempdir(), "filtdir"),
            rrna_db_dir = dbdir
        )
    }
    
    expect_equal(class(r1), "data.frame")
})

