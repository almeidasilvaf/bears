
#' Index the transcriptome for kallisto
#' 
#' @param kallistoindex Directory where kallisto index file will be stored.
#' Default: results/05_quantification/kallisto/idx.
#' @param transcriptome_path Path to the reference transcriptome FASTA file.
#'
#' @return A 2-column data frame with path to index file in the first column
#' and index build status in the second column, with "OK" if the transcriptome
#' index was successfully created, and NA otherwise.
#' @export
#' @rdname kallisto_index
#' @examples
#' kallistoindex <- file.path(tempdir(), "idx")
#' transcriptome_path <- system.file(
#'      "extdata", "Hsapiens_GRCh37.75_subset_transcripts.fa", package="bears"
#' )
#' if(kallisto_is_installed()) {
#'     kallisto_index(kallistoindex, transcriptome_path)
#' }
kallisto_index <- function(
    kallistoindex = "results/05_quantification/kallisto/idx",
    transcriptome_path = NULL) {
    
    if(!kallisto_is_installed()) { stop("Unable to find kallisto in PATH.") }
    if(!dir.exists(kallistoindex)) { 
        dir.create(kallistoindex, recursive = TRUE) 
    }
    idx <- paste0(kallistoindex, "/transcripts.idx")
    args <- c("index -i", idx, transcriptome_path)
    system2("kallisto", args = args)
    status <- "OK"
    if(!file.exists(idx)) { status <- NA }
    df_status <- data.frame(index_path = idx, status = status)
    return(df_status)
}


#' Wrapper to handle technical replicates during kallisto quantification
#' 
#' @param sample_info Data frame of sample metadata created with the
#' functions \code{create_sample_info} and \code{infer_strandedness}.
#' @param filtdir Path to the directory where filtered reads are stored.
#' Default: results/03_filtered_FASTQ.
#' 
#' @return A list with 2 elements: 
#' \describe{
#'   \item{paired}{List of BioSamples and paths to read pairs in a character
#'   object.}
#'   \item{single}{List of BioSamples and paths to reads in a character object.}
#' }
#' @noRd
run2biosample_kallisto <- function(sample_info = NULL, 
                                   filtdir = "results/03_filtered_FASTQ") {
    single <- sample_info[sample_info$Layout == "SINGLE", ]
    paired <- sample_info[sample_info$Layout == "PAIRED", ]
    
    if(nrow(paired) > 0) {
        paired$Run1 <- paste0(filtdir, "/", paired$Run, "_1.fastq.gz")
        paired$Run2 <- paste0(filtdir, "/", paired$Run, "_2.fastq.gz")
        pair_list <- split(paired, paired$BioSample)
        pair_list <- lapply(pair_list, function(x) {
            x$pair <- paste(x$Run1, x$Run2, sep = " ")
            pairs <- paste(x$pair, collapse = " ")
            return(pairs)
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


#' Quantify expression with kallisto
#' 
#' @param sample_info Data frame of sample metadata created with the
#' functions \code{create_sample_info} and \code{infer_strandedness}.
#' The function \code{infer_strandedness} adds a column named "Orientation" 
#' with library strandedness information, which is mandatory for 
#' kallisto quantification.
#' @param qc_table Data frame of fastp summary statistics as returned
#' by \code{summary_stats_fastp()}.
#' @param filtdir Path to the directory where filtered reads are stored.
#' Default: results/03_filtered_FASTQ.
#' @param kallistoindex Directory where kallisto index file will be stored.
#' Default: results/05_quantification/kallisto/idx.
#' @param kallistodir Directory where quantification files will be stored.
#' Default: results/05_quantification/kallisto.
#' @param threads Number of threads for kallisto quant.
#'
#' @return A 2-column data frame with BioSample IDs in the first column and 
#' quantification status in the second column, with "OK" if kallisto 
#' successfully quantified expression for a given BioSample, and NA otherwise.
#' @export
#' @rdname kallisto_quantify
#' @examples
#' data(sample_info)
#' qc_table <- summary_stats_fastp(system.file("extdata", package = "bears"))
#' 
#' filtdir <- system.file("extdata", package = "bears")
#' kallistoindex <- file.path(tempdir(), "transcripts.idx")
#' kallistodir <- tempdir()
#' transcriptome_path <- system.file(
#'      "extdata", "Hsapiens_GRCh37.75_subset_transcripts.fa", package="bears"
#' )
#' if(kallisto_is_installed()) {
#'     kallisto_index(kallistoindex, transcriptome_path)
#'     kallisto_quantify(
#'         sample_info, qc_table, filtdir, kallistoindex, kallistodir
#'     )
#' }
kallisto_quantify <- function(
    sample_info = NULL, 
    qc_table = NULL,
    filtdir = "results/03_filtered_FASTQ",
    kallistoindex = "results/05_quantification/kallisto/idx",
    kallistodir = "results/05_quantification/kallisto",
    threads = NULL
    ) {
    if(!kallisto_is_installed()) { stop("Unable to find kallisto in PATH.") }
    if(!dir.exists(kallistodir)) { dir.create(kallistodir, recursive = TRUE) }
    
    idx <- paste0(kallistoindex, "/transcripts.idx")
    r <- run2biosample_kallisto(sample_info, filtdir)
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
            outdir <- paste0(kallistodir, "/", var$biosample)
            orientation <- sample_meta[x, "Orientation"]
            lib <- translate_strandedness(orientation, var$layout)$kallisto
            args <- c("quant", "-i", idx, "-o", outdir, lib, 
                      "--bootstrap-samples=100 --bias")
            if(!is.null(threads)) { 
                args <- c(args, paste0("--threads=", threads)) 
            }
            if(var$layout == "SINGLE") {
                frag_len <- 100
                read_len <- qc_table[qc_table$Sample == var$run,
                                         "after_meanlength"]
                if(read_len > 100) {
                    frag_len <- read_len + 100
                }
                reads <- r$single[[var$biosample]]
                args <- c(args, "--single -l", frag_len, "-s 20", reads)
            } else if(var$layout == "PAIRED") {
                reads <- r$paired[[var$biosample]]
                args <- c(args, reads)
            } else {
                message("Layout information not available.")
            }
            system2("kallisto", args = args)
        }
    })
    dirlist <- list.dirs(kallistodir, recursive = FALSE, full.names = FALSE)
    df_status <- data.frame(sample = sample_meta$BioSample)
    df_status$status <- ifelse(df_status$sample %in% dirlist, "OK", NA)
    return(df_status)
}


#' Create a SummarizedExperiment object from kallisto output
#' 
#' @param sample_info Data frame of sample metadata created with the
#' functions \code{create_sample_info}
#' @param level Character indicating to which level expression must be 
#' quantified in the SE object. One of "gene" (default), "transcript", 
#' or "both". For "both", the SE object will have two assays named "transcript"
#' and "gene".
#' @param kallistodir Directory where quantification files will be stored.
#' Default: results/05_quantification/kallisto.
#' @param tx2gene Data frame of correspondence between genes and transcripts, 
#' with gene IDs in the first column and transcript IDs in the second column.
#' Only required if level = 'gene' or 'both'. 
#'
#' @return A SummarizedExperiment object with gene/transcript expression
#' levels and sample metadata.
#' @importFrom tximport tximport summarizeToGene
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @rdname kallisto2se
#' @examples 
#' data(sample_info)
#' data(tx2gene)
#' kallistodir <- system.file("extdata", package="bears")
#' se_gene <- kallisto2se(sample_info, kallistodir = kallistodir, 
#'                        tx2gene = tx2gene)
kallisto2se <- function(sample_info = NULL, level="gene", 
                        kallistodir = "results/05_quantification/kallisto", 
                        tx2gene = NULL) {
    sample_meta <- sample_info[!duplicated(sample_info$BioSample), ]
    files <- file.path(kallistodir, sample_meta$BioSample, "abundance.tsv")
    names(files) <- paste0(sample_meta$BioSample)
    coldata <- data.frame(row.names = sample_meta$BioSample)
    coldata <- cbind(coldata, sample_meta[, !names(sample_meta) == "BioSample"])
    
    if(level == "gene") {
        exp <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene)
        final <- SummarizedExperiment::SummarizedExperiment(
            assays = list(gene_TPM = exp$abundance, gene_counts = exp$counts),
            colData = coldata
        )
    } else if(level == "transcript") {
        exp <- tximport::tximport(files, type = "kallisto", txOut = TRUE)
        final <- SummarizedExperiment::SummarizedExperiment(
            assays = list(tx_TPM = exp$abundance, tx_counts = exp$counts),
            colData = coldata
        )
    } else if(level == "both") {
        exp_tx <- tximport::tximport(files, type = "kallisto", txOut = TRUE)
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
