
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

