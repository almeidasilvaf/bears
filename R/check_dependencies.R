
#' Check prerequisites to load conda environment temporarily with Herper
#' 
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#'
#' @noRd
#' @return Logical indicating whether conda environment should be loaded on not.
load_env <- function(envname = NULL, miniconda_path = NULL) {
    if(!is.null(envname) & !is.null(miniconda_path)) {
        load <- TRUE
    } else if(!is.null(envname) & is.null(miniconda_path)) {
        stop("To load a conda environment, both `envname` and `miniconda_path`
             must be defined.")
    } else if(is.null(envname) & !is.null(miniconda_path)) {
        stop("To load a conda environment, both `envname` and `miniconda_path`
             must be defined.")
    } else {
        load <- FALSE
    }
    return(load)
}


#' Wrapper to check if command is found in PATH
#' 
#' @param cmd Command to test.
#' @param args Arguments for command.
#'
#' @return Logical indicating whether the command is in PATH or not.
#' @noRd
is_valid <- function(cmd = NULL, args = NULL) {
    found <- tryCatch(
        system2(cmd, args = args, stdout = FALSE, stderr = FALSE), 
        error = function(e) return(FALSE),
        warning = function(w) return(FALSE)
    )
    if(!isFALSE(found)) {
        found <- TRUE
    }
    return(found)
}

#' Check if SRAToolkit is installed
#' 
#' @return Logical indicating whether SRAToolkit is installed or not.
#' @export
#' @rdname sratoolkit_is_installed
#' @examples 
#' sratoolkit_is_installed()
sratoolkit_is_installed <- function() {
    valid <- is_valid(cmd = "fastq-dump")
    return(valid)
}


#' Check if FastQC is installed
#' 
#' @return Logical indicating whether FastQC is installed or not.
#' @export
#' @rdname fastqc_is_installed
#' @examples 
#' sratoolkit_is_installed()
fastqc_is_installed <- function() {
    valid <- is_valid(cmd = "fastqc", args = "-h")
    return(valid)
}


#' Check if multiqc is installed
#' 
#' @return Logical indicating whether MultiQC is installed or not.
#' @export
#' @rdname multiqc_is_installed
#' @examples 
#' multiqc_is_installed()
multiqc_is_installed <- function() {
    valid <- is_valid(cmd = "multiqc", args = "-h")
    return(valid)
}


#' Check if Trimmomatic is installed
#' 
#' @return Logical indicating whether Trimmomatic is installed or not.
#' @export
#' @rdname trimmomatic_is_installed
#' @examples 
#' trimmomatic_is_installed()
trimmomatic_is_installed <- function() {
    valid <- is_valid(cmd = "trimmomatic", args = "-version")
    return(valid)
}


#' Check if SortMeRNA is installed
#' 
#' @return Logical indicating whether SortMeRNA is installed or not.
#' @export
#' @rdname sortmerna_is_installed
#' @examples 
#' sortmerna_is_installed()
sortmerna_is_installed <- function() {
    valid <- is_valid(cmd = "sortmerna", args = "-h")
    return(valid)
}


#' Check if STAR is installed
#' 
#' @return Logical indicating whether STAR is installed or not.
#' @export
#' @rdname star_is_installed
#' @examples
#' star_is_installed()
star_is_installed <- function() {
    valid <- is_valid(cmd = "STAR", args = "-h")
    return(valid)
}


#' Check if SHRiMP is installed
#' 
#' @return Logical indicating whether SHRiMP is instaled or not
#' @export
#' @rdname shrimp_is_installed
#' @examples 
#' shrimp_is_installed()
shrimp_is_installed <- function() {
    valid <- is_valid(cmd = "gmapper-cs", args = "-h")
    return(valid)
}


#' Check if RSeQC is installed
#' 
#' @return Logical indicating whether RSeQC is instaled or not
#' @export
#' @rdname rseqc_is_installed
#' @examples 
#' rseqc_is_installed()
rseqc_is_installed <- function() {
    valid <- is_valid(cmd = "infer_experiment.py", args = "-h")
    return(valid)
}


#' Check if salmon is installed
#' 
#' @return Logical indicating whether salmon is instaled or not
#' @export
#' @rdname salmon_is_installed
#' @examples 
#' salmon_is_installed()
salmon_is_installed <- function() {
    valid <- is_valid(cmd = "salmon", args = "-h")
    return(valid)
}









