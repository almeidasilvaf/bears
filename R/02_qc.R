
#' Wrapper to get paths to .fastq files depending on the layout
#' 
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' @param run Run accession.
#' @param cmd What command to run. One of "S", "P1", or "P2". S will
#' run the command for single-end reads. P1 will run the command for 
#' the first pair of paired-end reads. P2 will run the command for 
#' the second pair.
#' @return Character with the path to .fastq file.
#' @noRd
get_fastq_paths <- function(fastqdir, run, cmd = "S") {
    fastqfile <- NULL
    if(cmd == "S") {
        if(file.exists(paste0(fastqdir, "/", run, ".fastq.gz"))) {
            fastqfile <- paste0(fastqdir, "/", run, ".fastq.gz")
        } else { message("File not found") }
    } else if(cmd == "P1") {
        if(file.exists(paste0(fastqdir, "/", run, "_1.fastq.gz"))) {
            fastqfile <- paste0(fastqdir, "/", run, "_1.fastq.gz")
        } else if(file.exists(paste0(fastqdir, "/", run, ".1.fastq.gz"))) {
            fastqfile <- paste0(fastqdir, "/", run, ".1.fastq.gz")
        } else { message("File not found.") }
    } else {
        if(file.exists(paste0(fastqdir, "/", run, "_2.fastq.gz"))) {
            fastqfile <- paste0(fastqdir, "/", run, "_2.fastq.gz")
        } else if(file.exists(paste0(fastqdir, "/", run, ".2.fastq.gz"))) {
            fastqfile <- paste0(fastqdir, "/", run, ".2.fastq.gz")
        } else { message("File not found.") }
    }
    return(fastqfile)
}


#' Wrapper to get paths to run FastQC depending on the layout
#' 
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' @param fastqcdir Path to the directory where FastQC output files will
#' be stored.
#' @param run Run accession.
#' @param cmd What command to run. One of "S", "P1", or "P2". S will
#' run the command for single-end reads. P1 will run the command for 
#' the first pair of paired-end reads. P2 will run the command for 
#' the second pair.
#'  
#' @return Paths to fastqfile, old HTML, new HTML, old .zip, and new .zip.
#' @noRd
get_fastqc_paths <- function(fastqdir = NULL, fastqcdir = NULL,
                             run = NULL, cmd = "S") {
    fastqfile <- get_fastq_paths(fastqdir, run, cmd)
    if(cmd == "S") {
        if(file.exists(paste0(fastqdir, "/", run, ".fastq.gz"))) {
            old_htmlpath <- paste0(fastqdir, "/", run, "_fastqc.html")
            new_htmlpath <- paste0(fastqcdir, "/", run, "_fastqc.html")
            old_zippath <- paste0(fastqdir, "/", run, "_fastqc.zip")
            new_zippath <- paste0(fastqcdir, "/", run, "_fastqc.zip")
        } else { message("File not found") }
    } else if(cmd == "P1") {
        if(file.exists(paste0(fastqdir, "/", run, "_1.fastq.gz"))) {
            old_htmlpath <- paste0(fastqdir, "/", run, "_1_fastqc.html")
            new_htmlpath <- paste0(fastqcdir, "/", run, "_1_fastqc.html")
            old_zippath <- paste0(fastqdir, "/", run, "_1_fastqc.zip")
            new_zippath <- paste0(fastqcdir, "/", run, "_1_fastqc.zip")
        } else if(file.exists(paste0(fastqdir, "/", run, ".1.fastq.gz"))) {
            old_htmlpath <- paste0(fastqdir, "/", run, ".1_fastqc.html")
            new_htmlpath <- paste0(fastqcdir, "/", run, ".1_fastqc.html")
            old_zippath <- paste0(fastqdir, "/", run, ".1_fastqc.zip")
            new_zippath <- paste0(fastqcdir, "/", run, ".1_fastqc.zip")
        } else { message("File not found.") }
    } else {
        if(file.exists(paste0(fastqdir, "/", run, "_2.fastq.gz"))) {
            old_htmlpath <- paste0(fastqdir, "/", run, "_2_fastqc.html")
            new_htmlpath <- paste0(fastqcdir, "/", run, "_2_fastqc.html")
            old_zippath <- paste0(fastqdir, "/", run, "_2_fastqc.zip")
            new_zippath <- paste0(fastqcdir, "/", run, "_2_fastqc.zip")
        } else if(file.exists(paste0(fastqdir, "/", run, ".2.fastq.gz"))) {
            old_htmlpath <- paste0(fastqdir, "/", run, ".2_fastqc.html")
            new_htmlpath <- paste0(fastqcdir, "/", run, ".2_fastqc.html")
            old_zippath <- paste0(fastqdir, "/", run, ".2_fastqc.zip")
            new_zippath <- paste0(fastqcdir, "/", run, ".2_fastqc.zip")
        } else { message("File not found.") }
    }
    path <- list(fastq = fastqfile, 
                 oldhtml = old_htmlpath, newhtml = new_htmlpath,
                 oldzip = old_zippath, newzip = new_zippath)
    return(path)
}

#' Run FastQC on .fastq files
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' Default: results/01_FASTQ_files.
#' @param fastqcdir Path to the directory where FastQC output will be stored.
#' Default: results/02_FastQC_dir.
#' 
#' @importFrom fs file_move
#' @return A 2-column data frame with run accession in the first column and
#' FastQC run status. If FastQC ran successfully, the status "OK" is displayed.
#' If FastQC failed to process a particular file, NA will be displayed on its 
#' corresponding name.
#' @export
#' @rdname run_fastqc
#' @examples
#' data(sample_info)
#' fq <- system.file("extdata", package="bears")
#' fqc <- tempdir()
#' if(fastqc_is_installed()) {
#'     run_fastqc(sample_info, fastqdir = fq, fastqcdir = fqc)
#' }
run_fastqc <- function(sample_info,
                       fastqdir = "results/01_FASTQ_files",
                       fastqcdir = "results/02_FastQC_dir") {
    if(!dir.exists(fastqcdir)) { dir.create(fastqcdir, recursive = TRUE) }
    if(!fastqc_is_installed()) { stop("Unable to find FastQC in PATH.") }
    
    d <- lapply(seq_len(nrow(sample_info)), function(x) {
        var <- var2list(sample_info, index = x)
        file <- paste0(fastqdir, "/", var$run, ".fastq.gz")
        if(var$layout == "PAIRED") { 
            file <- paste0(fastqdir, "/", var$run, "_1.fastq.gz") 
        }
        if(skip(var$platform, path = file)) {
            message("Skipping file...")
        } else {
            if(var$layout == "SINGLE") {
                p <- get_fastqc_paths(fastqdir, fastqcdir, var$run, cmd="S")
                system2("fastqc", args = p$fastq)
                move_html <- fs::file_move(p$oldhtml, p$newhtml)
                move_zip <- fs::file_move(p$oldzip, p$newzip)
            } else if(var$layout == "PAIRED") {
                p1 <- get_fastqc_paths(fastqdir, fastqcdir, var$run, cmd="P1")
                p2 <- get_fastqc_paths(fastqdir, fastqcdir, var$run, cmd="P2")
                system2("fastqc", args = p1$fastq)
                system2("fastqc", args = p2$fastq)
                move_html1 <- fs::file_move(p1$oldhtml, p1$newhtml)
                move_html2 <- fs::file_move(p2$oldhtml, p2$newhtml)
                move_zip1 <- fs::file_move(p1$oldzip, p1$newzip)
                move_zip2 <- fs::file_move(p2$oldzip, p2$newzip)
            } else {
                message("Layout information not available.")
            }
        }
    })
    flist <- list.files(path = fastqcdir, pattern = ".zip")
    flist <- unique(gsub("_[0-9]_fastqc.zip|_fastqc.zip", "", flist))
    df_status <- data.frame(run = sample_info$Run)
    df_status$status <- ifelse(flist %in% df_status$run, "OK", NA)
    return(df_status)
}

#' Run MultiQC to get a summary of QC results
#' 
#' @param dir Directory where MultiQC should be run. 
#' Default: results/02_FastQC_dir (for FastQC results).
#' @param outdir Path the output directory. 
#' Default: results/multiqc/fastqc.
#' @param runon Type of QC report to generate. One of "fastqc" or "star".
#' 
#' @return A data frame with QC summary for each sample.
#' @export
#' @rdname multiqc
#' @examples
#' data(sample_info)
#' fq <- system.file("extdata", package="bears")
#' out <- tempdir()
#' if(multiqc_is_installed()) {
#'     fastqc_table <- multiqc(fq, out)
#' }
multiqc <- function(dir="results/02_FastQC_dir",
                    outdir = "results/multiqc/fastqc",
                    runon="fastqc") {
    if(!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
    if(!multiqc_is_installed()) { stop("Unable to find MultiQC in PATH.") }
    
    args <- c("-o", outdir, dir)
    system2("multiqc", args = args)
    if(runon == "fastqc") {
        fpath <- paste0(outdir, "/multiqc_data/multiqc_fastqc.txt")
    } else {
        fpath <- paste0(outdir, "/multiqc_data/multiqc_star.txt")
    }
    qc_table <- read.csv(fpath, header=TRUE, sep="\t")
    return(qc_table)
}







