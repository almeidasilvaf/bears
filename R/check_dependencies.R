

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


#' Check if FastQC is installed
#' 
#' @return Logical indicating whether FastQC is installed or not.
#' @export
#' @rdname fastqc_is_installed
#' @examples 
#' fastqc_is_installed()
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

#' Check if RSeQC is installed
#' 
#' @return Logical indicating whether RSeQC is installed or not
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
#' @return Logical indicating whether salmon is installed or not
#' @export
#' @rdname salmon_is_installed
#' @examples 
#' salmon_is_installed()
salmon_is_installed <- function() {
    valid <- is_valid(cmd = "salmon", args = "-h")
    return(valid)
}


#' Check if kallisto is installed
#' 
#' @return Logical indicating whether kallisto is installed or not
#' @export
#' @rdname kallisto_is_installed
#' @examples 
#' kallisto_is_installed()
kallisto_is_installed <- function() {
    valid <- is_valid(cmd = "kallisto", args = "-h")
    return(valid)
}


#' Check if subread is installed
#' 
#' @return Logical indicating whether subread is installed or not
#' @export
#' @rdname subread_is_installed
#' @examples 
#' subread_is_installed()
subread_is_installed <- function() {
    valid <- is_valid(cmd = "featureCounts", args = "-h")
    return(valid)
}


#' Check if StringTie is installed
#' 
#' @return Logical indicating whether StringTie is installed or not
#' @export
#' @rdname stringtie_is_installed
#' @examples
#' stringtie_is_installed()
stringtie_is_installed <- function() {
    valid <- is_valid(cmd = "stringtie", args = "-h")
    return(valid)
}


#' Check if TACO is installed
#' 
#' @return Logical indicating whether TACO is installed or not
#' @export
#' @rdname taco_is_installed
#' @examples
#' taco_is_installed()
taco_is_installed <- function() {
    valid <- is_valid(cmd = "taco_run", args = "-h")
    return(valid)
}






