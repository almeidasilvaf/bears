
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
    
    expect_false(s1)
    expect_true(s2)
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
