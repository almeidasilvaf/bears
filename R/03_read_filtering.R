

#' Trim low-quality bases and adapters using fastp
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' Default: results/01_FASTQ_files.
#' @param qcdir Character indicating the path to the directory where
#' output summary stats will be saved. Default: results/QC_dir/fastp_stats.
#' @param filtdir Path to the directory where filtered .fastq files will
#' be temporarily stored. After trimming, filtered reads are moved back to
#' fastqdir. Default: results/03_filtered_FASTQ.
#' @param threads Numeric indicating the number of threads to use in fastp.
#' Default: 1.
#' @param delete_raw Logical indicating whether to delete raw (unfiltered)
#' FASTQ files after QC and filtering with fastp. It is recommended to
#' delete raw files, as they use much memory and the useful information
#' will be on filtered files, even if no filtering is performed. 
#' Default: FALSE.
#' 
#' @export
#' @rdname trim_reads
#' @importFrom fs file_move file_delete
#' @return A 2-column data frame with run accession in the first column
#' and fastp run status in the second column, with "OK" if it ran 
#' successfully and NA if a file could not be run.
#' @examples
#' data(sample_info)
#' fastqdir <- tempdir()
#' file.copy(system.file("extdata", "SRR1039508_1.fastq.gz", package = "bears"),
#'           fastqdir)
#' file.copy(system.file("extdata", "SRR1039508_2.fastq.gz", package = "bears"),
#'           fastqdir)
#' filtdir <- paste0(fastqdir, "/filtdir")
#' qcdir <- file.path(tempdir(), "qcdir")
#' if(!dir.exists(filtdir)) { dir.create(filtdir, recursive = TRUE) }
#' if(!dir.exists(qcdir)) { dir.create(qcdir, recursive = TRUE) }
#' if(fastp_is_installed()) {
#'     trim_status <- trim_reads(sample_info, fastqdir, filtdir, qcdir)
#' }
trim_reads <- function(sample_info = NULL, 
                       fastqdir = "results/01_FASTQ_files",
                       filtdir = "results/03_filtered_FASTQ",
                       qcdir = "results/QC_dir/fastp_stats",
                       threads = 1,
                       delete_raw = FALSE) {
    
    if(!dir.exists(filtdir)) { dir.create(filtdir, recursive = TRUE) }
    if(!dir.exists(qcdir)) { dir.create(qcdir, recursive = TRUE) }
    if(!fastp_is_installed()) { stop("Unable to find fastp in PATH.") }
    
    if(nrow(sample_info) > 0) {
        t <- lapply(seq_len(nrow(sample_info)), function(x) {
            var <- var2list(sample_info, index = x)
            file <- paste0(fastqdir, "/", var$run, ".fastq.gz")
            if(var$layout == "PAIRED") { 
                file <- paste0(fastqdir, "/", var$run, "_1.fastq.gz") 
            }
            if(skip(var$platform, path = file)) {
                message("Skipping file...")
            } else {
                
                json <- file.path(qcdir, paste0(var$run, ".json"))
                html <- file.path(qcdir, paste0(var$run, ".html"))
                
                if(var$layout == "SINGLE") {
                    file <- get_fastq_paths(fastqdir, var$run, cmd = "S")
                    filtfile <- paste0(filtdir, "/", var$run, ".fastq.gz")
                    args <- c(
                        "--in1", file, "--out1", filtfile, "--thread", threads,
                        "--json", json, "--html", html
                    )
                    system2("fastp", args = args)
                    if(delete_raw) { dfile <- fs::file_delete(file) } 
                    
                } else {
                    p1 <- get_fastq_paths(fastqdir, var$run, cmd = "P1")
                    p2 <- get_fastq_paths(fastqdir, var$run, cmd = "P2")
                    fp1 <- paste0(filtdir, "/", var$run, "_1.fastq.gz")
                    fp2 <- paste0(filtdir, "/", var$run, "_2.fastq.gz")
                    args <- c(
                        "--in1", p1, "--in2", p2, "--thread", threads,
                        "--out1", fp1, "--out2", fp2, 
                        "--json", json, "--html", html,
                        "--detect_adapter_for_pe"
                    )
                    system2("fastp", args = args)
                    if(delete_raw) {
                        dfile1 <- fs::file_delete(p1)
                        dfile2 <- fs::file_delete(p2)
                    }
                } 
            }
        })
    }
    flist <- basename(list.files(path = qcdir, pattern = ".json$"))
    flist <- unique(gsub("\\.json", "", flist))
    df_status <- data.frame(run = sample_info$Run)
    df_status$status <- ifelse(flist %in% df_status$run, "OK", NA)
    return(df_status)
}

#' Wrapper to create list of arguments to handle rRNA dbs in SortMeRNA 
#' 
#' @param rrna_db_dir Path to directory containing reference rRNA database,
#' which must be stored as FASTA files.
#' 
#' @return Character vector with arguments to pass to SortMeRNA in 
#' the system2 call inside \code{remove_rrna()}.
#' @noRd
rrna_ref_args <- function(rrna_db_dir = NULL) {
    fasta_files <- list.files(rrna_db_dir, full.names = TRUE,
                              pattern = ".fasta$|.fa$")
    
    refs <- unlist(lapply(fasta_files, function(x) {
        arg <- c("--ref", x)
        return(arg)
    }))
    return(refs)
}


#' Wrapper to clean directory after running SortMeRNA
#' 
#' @param filtdir Path to the directory where filtered reads will be stored.
#' Default: results/03_filtered_FASTQ.
#' 
#' @return A NULL object.
#' @noRd
#' @importFrom fs file_move file_delete
clean_sortmerna <- function(filtdir = "results/03_filtered_FASTQ") {
    # Delete rRNA
    old_files <- list.files(filtdir, pattern = "*.gz", full.names=TRUE)
    rrna <- old_files[grep("rRNA", old_files)]
    old_files <- old_files[!(old_files %in% rrna)]
    del <- fs::file_delete(rrna)
    
    # Rename files
    new_files <- gsub("_filt_fwd.fq.gz", "_1.fastq.gz", old_files)
    new_files <- gsub("_filt_rev.fq.gz", "_2.fastq.gz", new_files)
    new_files <- gsub("_filt.fq.gz", ".fastq.gz", new_files)
    rename <- fs::file_move(old_files, new_files)
    return(NULL)
}


#' Remove rRNA sequences from .fastq files with SortMeRNA
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param fastqdir Path to the directory where .fastq files are stored.
#' @param filtdir Path to the directory where filtered reads will be stored.
#' Default: results/03_filtered_FASTQ.
#' @param rrna_db_dir Path to directory containing reference rRNA database,
#' which must be stored as FASTA files.
#' @param threads Number of threads for SortMeRna. Default: 1.
#' 
#' @return A 2-column data frame with run accessions in the first column
#' and SortMeRNA running status in the second column, with "OK" if SortMeRNA
#' ran successfully for each file and NA otherwise.
#' @export
#' @rdname remove_rrna
#' @importFrom fs dir_delete
#' @examples 
#' data(sample_info)
#' fastqdir <- system.file("extdata", package="bears")
#' filtdir <- tempdir()
#' rrna_db_dir <- tempdir()
#' rrna_file <- system.file("extdata", "bac_16s_subset.fa", package="bears")
#' file.copy(from = rrna_file, to = rrna_db_dir)
#' if(sortmerna_is_installed()) {
#'     remove_rrna(sample_info, fastqdir, filtdir, rrna_db_dir)
#' }
remove_rrna <- function(sample_info,
                        fastqdir = "results/01_FASTQ_files",
                        filtdir = "results/03_filtered_FASTQ",
                        rrna_db_dir = NULL,
                        threads = 1) {
    if(!dir.exists(filtdir)) { dir.create(filtdir, recursive = TRUE) }
    if(!sortmerna_is_installed()) { stop("Unable to find SortMeRNA in PATH.") }
    refs <- rrna_ref_args(rrna_db_dir)
    t <- lapply(seq_len(nrow(sample_info)), function(x) {
        var <- var2list(sample_info, index = x)
        file <- paste0(fastqdir, "/", var$run, ".fastq.gz")
        if(var$layout == "PAIRED") { 
            file <- paste0(fastqdir, "/", var$run, "_1.fastq.gz") 
        }
        if(skip(var$platform, path = file)) {
            message("Skipping file...")
        } else {
            workdir <- paste0(filtdir, "/", var$run)
            if(var$layout == "SINGLE") {
                r <- get_fastq_paths(fastqdir, var$run, cmd = "S")
                args <- c(refs, "--reads", r, "--threads", threads,
                          "--workdir", workdir, "--fastx",
                          "--aligned", paste0(workdir, "_rRNA"),
                          "--other", paste0(workdir, "_filt"))
                system2("sortmerna", args = args)
            } else {
                r1 <- get_fastq_paths(fastqdir, var$run, cmd = "P1")
                r2 <- get_fastq_paths(fastqdir, var$run, cmd = "P2")
                args <- c(refs, "--reads", r1, "--reads", r2,
                          "--threads", threads, 
                          "--workdir", workdir, "--fastx",
                          "--aligned", paste0(workdir, "_rRNA"),
                          "--other", paste0(workdir, "_filt"),
                          "--paired_in", "--out2")
                system2("sortmerna", args = args)
            }
            delete_workdir <- fs::dir_delete(workdir)
        }
    })
    logdir <- paste0(filtdir, "/logs_sortmerna")
    if(!dir.exists(logdir)) { dir.create(logdir, recursive = TRUE) }
    logfiles <- list.files(filtdir, pattern = "*.log", full.names = TRUE)
    move_logs <- lapply(logfiles, function(x) fs::file_move(x, logdir))
    clean <- clean_sortmerna(filtdir)
    flist <- list.files(path = filtdir, pattern=".fastq.gz")
    flist <- unique(gsub("_[0-9].fastq.gz|.fastq.gz", "", flist))
    df_status <- data.frame(run = sample_info$Run)
    df_status$status <- ifelse(flist %in% df_status$run, "OK", NA)
    return(df_status)
}



