
#----Load data------------------------------------------------------------------
data(sample_info)

## Path to files and directories used in this set of tests
mappingdir <- system.file("extdata", package = "bears")
gff_path <- system.file(
    "extdata", "Homo_sapiens.GRCh37.75_subset.gtf", package = "bears"
)
fcountsdir <- file.path(tempdir(), "fcountsdir")


#----Start tests----------------------------------------------------------------
test_that("featureCounts() and featureCounts2se() work", {
    
    fc <- matrix(NA, nrow = 10, ncol = 10)
    fc_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts), colData = colData
    )
    if(subread_is_installed()) {
        fc <- featureCounts(sample_info, mappingdir, gff_path, fcountsdir)
        fc_se <- featureCounts2se(sample_info, fc)
    }
    
    expect_true(is.matrix(fc))
    expect_true("SummarizedExperiment" %in% class(fc_se))
})


