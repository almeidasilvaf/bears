
#' Index the transcriptome for kallisto
#' 
#' @param kallistoindex Directory where kallisto index file will be stored.
#' Default: results/05_quantification/kallisto/idx.
#' @param transcriptome_path Path to the reference transcriptome FASTA file.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#'
#' @return A NULL object.
#' @export
#' @rdname kallisto_index
#' @examples 
#' kallistoindex <- file.path(tempdir(), "transcripts.idx")
#' transcriptome_path <- system.file(
#'      "extdata", "Hsapiens_GRCh37.75_subset_transcripts.fa", package="bears"
#' )
#' if(kallisto_is_installed()) {
#'     kallisto_index(kallistoindex, transcriptome_path)
#' }
kallisto_index <- function(
    kallistoindex = "results/05_quantification/kallisto/idx",
    transcriptome_path = NULL,
    envname = NULL, miniconda_path = NULL) {
    
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!kallisto_is_installed()) { stop("Unable to find kallisto in PATH.") }
    if(!dir.exists(kallistoindex)) { 
        dir.create(kallistoindex, recursive = TRUE) 
    }
    idx <- paste0(kallistoindex, "/transcripts.idx")
    args <- c("index -i", idx, transcriptome_path)
    system2("kallisto", args = args)
    return(NULL)
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
#' @param fastqc_table Data frame of summary statistics for FastQC as returned
#' by \code{multiqc()}.
#' @param filtdir Path to the directory where filtered reads are stored.
#' Default: results/03_filtered_FASTQ.
#' @param kallistoindex Directory where kallisto index file will be stored.
#' Default: results/05_quantification/kallisto/idx.
#' @param kallistodir Directory where quantification files will be stored.
#' Default: results/05_quantification/kallisto.
#' @param threads Number of threads for kallisto quant.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#'
#' @return A NULL object.
#' @export
#' @rdname kallisto_quantify
#' @examples
#' data(sample_info)
#' data(fastqc_table)
#' filtdir <- system.file("extdata", package = "bears")
#' kallistoindex <- file.path(tempdir(), "transcripts.idx")
#' kallistodir <- tempdir()
#' transcriptome_path <- system.file(
#'      "extdata", "Hsapiens_GRCh37.75_subset_transcripts.fa", package="bears"
#' )
#' if(kallisto_is_installed()) {
#'     kallisto_index(kallistoindex, transcriptome_path)
#'     kallisto_quantify(sample_info, fastqc_table, filtdir, kallistoindex,
#'                       kallistodir)
#' }
kallisto_quantify <- function(
    sample_info = NULL, 
    fastqc_table = NULL,
    filtdir = "results/03_filtered_FASTQ",
    kallistoindex = "results/05_quantification/kallisto/idx",
    kallistodir = "results/05_quantification/kallisto",
    threads = NULL,
    envname = NULL, miniconda_path = NULL
    ) {
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!kallisto_is_installed()) { stop("Unable to find kallisto in PATH.") }
    if(!dir.exists(kallistodir)) { dir.create(kallistodir, recursive = TRUE) }
    
    idx <- paste0(kallistoindex, "/transcripts.idx")
    r <- run2biosample_kallisto(sample_info, filtdir)
    sample_meta <- sample_info[!duplicated(sample_info$BioSample), ]
    t <- lapply(seq_len(nrow(sample_meta)), function(x) {
        var <- var2list(sample_meta, index = x)
        if(grepl("SOLiD|PacBio", var$platform)) {
            message("Skipping PacBio/SOLiD reads...")
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
                read_len <- fastqc_table[fastqc_table$Sample == var$biosample,
                                         "Sequence.length"]
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
    return(NULL)
}

