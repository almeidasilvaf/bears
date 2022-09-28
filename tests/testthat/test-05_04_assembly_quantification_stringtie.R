
#----Load data------------------------------------------------------------------
data(sample_info)
data(tx2gene)

qc_table <- summary_stats_fastp(system.file("extdata", package = "bears"))
mappingdir <- system.file("extdata", package = "bears")
gff_path <- system.file(
    "extdata", "Homo_sapiens.GRCh37.75_subset.gtf", package = "bears"
)
stringtiedir <- file.path(tempdir(), "stringtiedir")

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
test_that("transcriptome assembly and quantification with StringTie works()", {
    
    sa <- data.frame()
    tm <- data.frame()
    stq <- data.frame()
    se_gene <- se_random
    se_tx <- se_random
    se_both <- list()
    
    
    if(stringtie_is_installed() & taco_is_installed()) {
        ## Assemble transcriptomes with StringTie
        sa <- stringtie_assemble(
            sample_info, qc_table, mappingdir, gff_path, stringtiedir
        )
        
        ## Merge transcriptome assemblies with TACO
        tm <- taco_merge(sample_info, stringtiedir)
        
        ## Quantify the expression with StringTie
        stq <- stringtie_quantify(
            sample_info, qc_table, mappingdir, gff_path, stringtiedir
        )
        
        ## Parse StringTie output as a SummarizedExperiment object
        se_gene <- stringtie2se(sample_info, stringtiedir, "gene", tx2gene)
        se_tx <- stringtie2se(sample_info, stringtiedir, "transcript", tx2gene)
        se_both <- stringtie2se(sample_info, stringtiedir, "both", tx2gene)
        
        expect_error(
            stringtie2se(sample_info, stringtiedir, "error", tx2gene)
        )
    }
    
    expect_equal(class(sa), "data.frame")
    expect_equal(class(tm), "data.frame")
    
    expect_true("SummarizedExperiment" %in% class(se_gene))
    expect_true("SummarizedExperiment" %in% class(se_tx))
    expect_equal(class(se_both), "list")
})