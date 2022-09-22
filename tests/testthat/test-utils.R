
#----Load data------------------------------------------------------------------


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

test_that("")