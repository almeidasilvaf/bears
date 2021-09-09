
#' Index the transcriptome for salmon
#' 
#' @param salmonindex Directory where the transcriptome index will be stored.
#' Default: results/05_quantification/salmon/idx.
#' @param transcriptome_path Path to the reference transcriptome FASTA file.
#' @param klen K-mer length. Default: 31.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#'
#' @return A 2-column data frame with path to index in the first column and
#' index build status in the second column, with "OK" if the transcriptome
#' index was successfully created, and NA otherwise.
#' @export
#' @rdname salmon_index
#' @examples
#' salmonindex <- tempdir()
#' transcriptome_path <- system.file(
#'      "extdata", "Hsapiens_GRCh37.75_subset_transcripts.fa", package="bears"
#' )
#' if(salmon_is_installed()) {
#'     salmon_index(salmonindex, transcriptome_path)
#' }
salmon_index <- function(salmonindex = "results/05_quantification/salmon/idx",
    transcriptome_path = NULL, 
    klen = 31,
    envname = NULL,
    miniconda_path = NULL) {
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!salmon_is_installed()) { stop("Unable to find salmon in PATH.") }
    if(!dir.exists(salmonindex)) { dir.create(salmonindex, recursive = TRUE) }
    args <- c("index -t", transcriptome_path, "-i", salmonindex, "-k", klen)
    system2("salmon", args = args)
    status <- "OK"
    if(length(dir(salmonindex)) == 0) { status <- NA }
    df_status <- data.frame(index_path = salmonindex, status = status)
    return(df_status)
}

#' Wrapper to handle technical replicates during salmon quantification
#' 
#' @param sample_info Data frame of sample metadata created with the
#' functions \code{create_sample_info} and \code{infer_strandedness}.
#' The column "Orientation", added by \code{infer_strandedness}, is mandatory.
#' @param filtdir Path to the directory where filtered reads are stored.
#' Default: results/03_filtered_FASTQ.
#' 
#' @return A list with 2 elements: 
#' \describe{
#'   \item{paired}{List of lists, with each element corresponding to a 
#'   BioSample. For each BioSample, there are elements 'run1' and 'run2',
#'   which contain the paths to .fastq files.}
#'   \item{single}{List of lists, with each element corresponding to a
#'   BioSample. For each BioSample, there is a character object with 
#'   path to each .fastq file}
#' }
#' @noRd
run2biosample_salmon <- function(sample_info = NULL, 
                                 filtdir = "results/03_filtered_FASTQ") {
    single <- sample_info[sample_info$Layout == "SINGLE", ]
    paired <- sample_info[sample_info$Layout == "PAIRED", ]
    
    if(nrow(paired) > 0) {
        paired$Run1 <- paste0(filtdir, "/", paired$Run, "_1.fastq.gz")
        paired$Run2 <- paste0(filtdir, "/", paired$Run, "_2.fastq.gz")
        pair_list <- split(paired, paired$BioSample)
        pair_list <- lapply(pair_list, function(x) {
            y <- list(run1 = paste(x$Run1, collapse = " "), 
                      run2 = paste(x$Run2, collapse = " "))
            return(y)
        })
    }
    if(nrow(single) > 0) {
        single$Run <- paste0(filtdir, "/", single$Run, ".fastq.gz")
        single_list <- split(single, single$BioSample)
        single_list <- lapply(single_list, function(x) {
            y <- paste(x$Run, collapse = " ")
            return(y)
        })
    }
    if(!exists("pair_list")) { pair_list <- NULL }
    if(!exists("single_list")) { single_list <- NULL }
    final_list <- list(paired = pair_list, single = single_list)
    return(final_list)
}


#' Quantify expression with salmon
#' 
#' @param sample_info Data frame of sample metadata created with the
#' functions \code{create_sample_info} and \code{infer_strandedness}.
#' The function \code{infer_strandedness} adds a column named "Orientation" 
#' with library strandedness information. If this column is not present in 
#' sample_info, salmon will automatically infer the library strandedness, but
#' it will take longer to run.
#' @param filtdir Path to the directory where filtered reads are stored.
#' Default: results/03_filtered_FASTQ.
#' @param salmonindex Directory where the transcriptome index is stored.
#' Default: results/05_quantification/salmon/idx.
#' @param salmondir Directory where quantification files will be stored.
#' Default: results/05_quantification/salmon.
#' @param threads Number of threads for salmon quant.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#'
#' @return A 2-column data frame with BioSample IDs in the first column and
#' salmon quantification status in the second column, with "OK" if salmon
#' sucessfully quantified expression for a given BioSample, and NA otherwise.
#' @export
#' @rdname salmon_quantify
#' @examples
#' data(sample_info)
#' filtdir <- system.file("extdata", package = "bears")
#' salmonindex <- tempdir()
#' salmondir <- tempdir()
#' transcriptome_path <- system.file(
#'      "extdata", "Hsapiens_GRCh37.75_subset_transcripts.fa", package="bears"
#' )
#' if(salmon_is_installed()) {
#'     salmon_index(salmonindex, transcriptome_path)
#'     salmon_quantify(sample_info, filtdir, salmonindex, salmondir)
#' }
salmon_quantify <- function(sample_info = NULL,
                            filtdir = "results/03_filtered_FASTQ",
                            salmonindex = "results/05_quantification/salmon/idx",
                            salmondir = "results/05_quantification/salmon",
                            threads = NULL,
                            envname = NULL,
                            miniconda_path = NULL) {
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!salmon_is_installed()) { stop("Unable to find salmon in PATH.") }
    if(!dir.exists(salmondir)) { dir.create(salmondir, recursive = TRUE) }
    
    r <- run2biosample_salmon(sample_info, filtdir)
    sample_meta <- sample_info[!duplicated(sample_info$BioSample), ]
    t <- lapply(seq_len(nrow(sample_meta)), function(x) {
        var <- var2list(sample_meta, index = x)
        file <- paste0(filtdir, "/", var$run, ".fastq.gz")
        if(var$layout == "PAIRED") { 
            file <- paste0(filtdir, "/", var$run, "_1.fastq.gz") 
        }
        if(skip(var$platform, path = file)) {
            message("Skipping file...")
        } else {
            if(var$layout == "SINGLE") {
                read_arg <- c("-r", r$single[[var$biosample]])
            } else if(var$layout == "PAIRED") {
                read_arg <- c("-1", r$paired[[var$biosample]]$run1,
                              "-2", r$paired[[var$biosample]]$run2)
            } else {
                message("Layout information not available.")
            }
            orientation <- sample_meta[x, "Orientation"] 
            if("Orientation" %in% names(sample_info)) {
                lib <- translate_strandedness(orientation, var$layout)$salmon
            } else {
                lib <- "A"
            }
            outdir <- paste0(salmondir, "/", var$biosample)
            args <- c("quant -i", salmonindex, "-l", lib, read_arg,
                      "-o", outdir, "--seqBias --gcBias --dumpEq")
            if(!is.null(threads)) { args <- c(args, "-p", threads) }
            system2("salmon", args = args)
        }
    })
    dirlist <- list.dirs(salmondir, recursive = FALSE, full.names = FALSE)
    df_status <- data.frame(sample = sample_meta$BioSample)
    df_status$status <- ifelse(df_status$sample %in% dirlist, "OK", NA)
    return(df_status)
}


#' Create a SummarizedExperiment object from salmon output
#' 
#' @param sample_info Data frame of sample metadata created with the
#' functions \code{create_sample_info}
#' @param level Character indicating to which level expression must be 
#' quantified in the SE object. One of "gene" (default), "transcript", 
#' or "both". For "both", the SE object will have two assays named "transcript"
#' and "gene".
#' @param salmondir Directory where quantification files will be stored.
#' Default: results/05_quantification/salmon.
#' @param tx2gene Data frame of correspondence between genes and transcripts, 
#' with gene IDs in the first column and transcript IDs in the second column.
#' Only required if level = 'gene' or 'both'. 
#'
#' @return A SummarizedExperiment object with gene/transcript expression
#' levels and sample metadata.
#' @importFrom tximport tximport summarizeToGene
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @rdname salmon2se
#' @examples 
#' data(sample_info)
#' data(tx2gene)
#' salmondir <- system.file("extdata", package="bears")
#' se_gene <- salmon2se(sample_info, salmondir = salmondir, tx2gene = tx2gene)
salmon2se <- function(sample_info = NULL, level="gene", 
                      salmondir = "results/05_quantification/salmon", 
                      tx2gene = NULL) {
    sample_meta <- sample_info[!duplicated(sample_info$BioSample), ]
    files <- file.path(salmondir, sample_meta$BioSample, "quant.sf")
    names(files) <- paste0(sample_meta$BioSample)
    coldata <- data.frame(row.names = sample_meta$BioSample)
    coldata <- cbind(coldata, sample_meta[, !names(sample_meta) == "BioSample"])
    
    if(level == "gene") {
        exp <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene)
        final <- SummarizedExperiment::SummarizedExperiment(
            assays = list(gene_TPM = exp$abundance, gene_counts = exp$counts),
                          colData = coldata
        )
    } else if(level == "transcript") {
        exp <- tximport::tximport(files, type = "salmon", txOut = TRUE)
        final <- SummarizedExperiment::SummarizedExperiment(
            assays = list(tx_TPM = exp$abundance, tx_counts = exp$counts),
            colData = coldata
        )
    } else if(level == "both") {
        exp_tx <- tximport::tximport(files, type = "salmon", txOut = TRUE)
        exp_gene <- tximport::summarizeToGene(exp_tx, tx2gene)
        se_gene <- SummarizedExperiment::SummarizedExperiment(
            assays = list(gene_TPM = exp_gene$abundance, 
                          gene_counts = exp_gene$counts),
            colData = coldata
        )
        se_tx <- SummarizedExperiment::SummarizedExperiment(
            assays = list(tx_TPM = exp_tx$abundance, 
                          tx_counts = exp_tx$counts),
            colData = coldata
        )
        final <- list(gene = se_gene, transcript = se_tx)
    } else {
        stop("Invalid parameter for the 'level' argument.")
    }
    return(final)
}












