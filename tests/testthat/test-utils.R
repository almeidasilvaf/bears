
#----Load data------------------------------------------------------------------
data(sample_info)
mapping_passed <- sample_info
bedpath <- system.file(
    "extdata", "Homo_sapiens.GRCh37.75_subset.bed", package = "bears"
)
mappingdir <- system.file("extdata", package = "bears")

#----Start tests----------------------------------------------------------------
test_that("skip() works", {
    s1 <- skip(platform = "Illumina Hiseq")
    s2 <- skip(platform = "PacBio")
    s3 <- skip(platform = "PacBio", path = bedpath)
    
    expect_false(s1)
    expect_true(s2)
    expect_equal(class(s3), "logical")
})

test_that("gff2bed() converts a GFF file to BED", {
    gff_path <- system.file(
        "extdata", "Homo_sapiens.GRCh37.75_subset.gtf", package="bears"
    )
    gffdir <- tempdir()
    file.copy(from = gff_path, to=gffdir)
    gff_file <- list.files(gffdir, full.names = TRUE, pattern = ".gtf")
    
    g1 <- gff2bed(gff_file)
    
    expect_equal(class(g1), "character")
})


test_that("infer_strandedness() infers strandedness", {
    
    s1 <- data.frame()
    if(rseqc_is_installed()) {
        s1 <- infer_strandedness(mapping_passed, bedpath, mappingdir)
    }
    
    expect_equal(class(s1), "data.frame")
})


test_that("translate_strandedness() correctly creates params for programs", {
    
    t1 <- translate_strandedness("first", layout = "SINGLE")
    t2 <- translate_strandedness("second", layout = "PAIRED")
    t3 <- translate_strandedness("unstranded", layout = "SINGLE")
    t4 <- translate_strandedness("random", layout = "SINGLE")
    
    expect_equal(class(t1), "list")
    expect_equal(class(t2), "list")
    expect_equal(class(t3), "list")
    expect_equal(class(t4), "list")
})


test_that("c_createdir() creates a dir if it doesn't exist", {
    
    dir <- file.path(tempdir(), paste0("testdir", sample(1:1000, 1)))
    c <- c_createdir(dir)
    
    expect_equal(class(c), "logical")
})


test_that("create_dir_structure() returns a list of paths", {
    
    rootdir <- tempdir()
    ds <- create_dir_structure(rootdir)
    
    expect_equal(class(ds), "list")
    expect_equal(class(ds$fastqdir), "character")
})




