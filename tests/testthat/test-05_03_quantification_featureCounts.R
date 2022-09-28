
#----Load data------------------------------------------------------------------
data(sample_info)

## Path to files and directories used in this set of tests
mappingdir <- system.file("extdata", package = "bears")
gff_path <- system.file(
    "extdata", "Homo_sapiens.GRCh37.75_subset.gtf", package = "bears"
)
fcountsdir <- file.path(tempdir(), "fcountsdir")

## Simulate a SummarizedExperiment object
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
test_that("featureCounts() and featureCounts2se() work", {
    
    fc <- matrix(NA, nrow = 10, ncol = 10)
    fc_se <- se_random
    
    if(subread_is_installed()) {
        fc <- featureCounts(sample_info, mappingdir, gff_path, fcountsdir)
        fc_se <- featureCounts2se(sample_info, fc)
    }
    
    expect_true(is.matrix(fc))
    expect_true("SummarizedExperiment" %in% class(fc_se))
})


