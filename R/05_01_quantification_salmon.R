


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
#' @return A NULL object.
#' @importFrom tools file_ext
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
    return(NULL)
}


#' Quantify expression with salmon
#' 
#' @param sample_info Data frame of sample metadata created with the
#' functions \code{create_sample_info} and \code{infer_strandedness}.
#' The column "Orientation", added by \code{infer_strandedness}, is mandatory.
#' @param fastqc_table Data frame of summary statistics for FastQC as returned
#' by \code{multiqc()}.
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
#' @return A NULL object.
#' @importFrom tools file_ext
#' @export
#' @rdname salmon_quantify
#' @examples
#' data(sample_info)
#' data(fastqc_table)
#' filtdir <- system.file("extdata", package = "bears")
#' salmonindex <- tempdir()
#' salmondir <- tempdir()
#' transcriptome_path <- system.file(
#'      "extdata", "Hsapiens_GRCh37.75_subset_transcripts.fa", package="bears"
#' )
#' if(salmon_is_installed()) {
#'     salmon_index(salmonindex, transcriptome_path)
#'     salmon_quantify(sample_info, fastqc_table, filtdir, 
#'                     salmonindex, salmondir)
#' }
#' 
salmon_quantify <- function(sample_info = NULL,
                            fastqc_table = NULL,
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
    
    t <- lapply(seq_len(nrow(sample_info)), function(x) {
        var <- var2list(sample_info, index = x)
        if(grepl("SOLiD|PacBio", var$platform)) {
            message("Skipping PacBio/SOLiD reads...")
        } else {
            if(var$layout == "SINGLE") {
                r <- get_fastq_paths(filtdir, var$run, cmd = "S")
                read_arg <- c("-r", r)
            } else if(var$layout == "PAIRED") {
                r1 <- get_fastq_paths(filtdir, var$run, cmd = "P1")
                r2 <- get_fastq_paths(filtdir, var$run, cmd = "P2")
                read_arg <- c("-1", r1, "-2", r2)
            } else {
                message("Layout information not available.")
            }
            orientation <- sample_info[x, "Orientation"] 
            libtype <- translate_strandedness(orientation, var$layout)$salmon
            outdir <- paste0(salmondir, "/", var$run)
            args <- c("quant -i", salmonindex, "-l", libtype, read_arg,
                      "-o", outdir, "--seqBias --gcBias --dumpEq")
            if(!is.null(threads)) { args <- c(args, "-p", threads) }
            system2("salmon", args = args)
        }
    })
    return(NULL)
}










