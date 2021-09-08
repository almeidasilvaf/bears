
#' Run Trimmomatic on reads that failed on FastQC check
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param fastqc_table Data frame of summary statistics for FastQC as returned
#' by \code{multiqc()}.
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' Default: results/01_FASTQ_files.
#' @param filtdir Path to the directory where filtered .fastq files will
#' be temporarily stored. After trimming, filtered reads are moved back to
#' fastqdir. Default: results/03_filtered_FASTQ.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#' 
#' @export
#' @rdname trim_reads
#' @importFrom fs file_move file_delete
#' @return A 2-column data frame with run accession in the first column
#' and Trimmomatic run status in the second column, with "OK" if it ran 
#' successfully and NA if a file could not be run.
#' @examples
#' data(sample_info)
#' data(fastqc_table)
#' fastqdir <- tempdir()
#' file.copy(system.file("extdata", "SRR1039508_1.fastq.gz", package = "bears"),
#'           fastqdir)
#' file.copy(system.file("extdata", "SRR1039508_2.fastq.gz", package = "bears"),
#'           fastqdir)
#' dir.create(paste0(fastqdir, "/filtdir"))
#' filtdir <- paste0(fastqdir, "/filtdir")
#' if(trimmomatic_is_installed()) {
#'     trim_reads(sample_info, fastqc_table, fastqdir, filtdir)
#' }
trim_reads <- function(sample_info = NULL,
                       fastqc_table = NULL,
                       fastqdir = "results/01_FASTQ_files",
                       filtdir = "results/03_filtered_FASTQ",
                       envname = NULL,
                       miniconda_path = NULL) {
    if(!dir.exists(filtdir)) { dir.create(filtdir, recursive = TRUE) }
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!trimmomatic_is_installed()) { stop("Unable to find Trimmomatic in PATH.") }
    pe_ad <- system.file("extdata", "PE_adapter.fa", package="bears")
    se_ad <- system.file("extdata", "SE_adapter.fa", package="bears")
    failed <- fastqc_table[fastqc_table$per_base_sequence_quality == "fail", 1]
    failed <- unique(gsub("_1$|_2$", "", failed))
    sample_info2 <- sample_info[sample_info$Run %in% failed, ]
    if(nrow(sample_info2) > 0) {
        t <- lapply(seq_len(nrow(sample_info2)), function(x) {
            var <- var2list(sample_info2, index = x)
            if(grepl("SOLiD|PacBio", var$platform)) {
                message("Skipping PacBio/SOLiD reads...")
            } else {
                trim <- c("LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
                if(var$layout == "SINGLE") {
                    file <- get_fastq_paths(fastqdir, var$run, cmd="S")
                    filtfile <- paste0(filtdir, "/", var$run, ".fastq.gz")
                    clip <- paste0("ILLUMINACLIP:", se_ad, ":2:30:10")
                    args <- c("SE -phred33", file, filtfile, clip, trim)
                    system2("trimmomatic", args = args)
                    move_file <- fs::file_move(filtfile, file)
                } else if(var$layout == "PAIRED") {
                    p1 <- get_fastq_paths(fastqdir, var$run, cmd="P1")
                    p2 <- get_fastq_paths(fastqdir, var$run, cmd="P2")
                    fp1 <- paste0(filtdir, "/", var$run, "_1.fastq.gz")
                    fp2 <- paste0(filtdir, "/", var$run, "_2.fastq.gz")
                    fp1_u <- paste0(filtdir, "/", var$run, "_1.unpaired.fq.gz")
                    fp2_u <- paste0(filtdir, "/", var$run, "_2.unpaired.fq.gz")
                    clip <- paste0("ILLUMINACLIP:", pe_ad, ":2:30:10")
                    args <- c("PE -phred33", p1, p2, fp1, fp1_u, fp2, fp2_u)
                    system2("trimmomatic", args = c(args, clip, trim))
                    move_file1 <- fs::file_move(fp1, p1)
                    move_file2 <- fs::file_move(fp2, p2)
                } else {
                    message("Layout information not available.")
                }
            }
        })
    }
    flist <- list.files(path = fastqdir, pattern=".fastq.gz")
    flist <- unique(gsub("_[0-9].fastq.gz|.fastq.gz", "", flist))
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
                              pattern = ".fasta|.fa")
    
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
#' @param threads Number of threads for SortMeRna. Default: 2.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
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
                        threads = 2,
                        envname = NULL,
                        miniconda_path = NULL) {
    if(!dir.exists(filtdir)) { dir.create(filtdir, recursive = TRUE) }
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!sortmerna_is_installed()) { stop("Unable to find SortMeRNA in PATH.") }
    refs <- rrna_ref_args(rrna_db_dir)
    t <- lapply(seq_len(nrow(sample_info)), function(x) {
        var <- var2list(sample_info, index = x)
        if(grepl("SOLiD|PacBio", var$platform)) {
            message("Skipping PacBio/SOLiD reads...")
        } else {
            workdir <- paste0(filtdir, "/", var$run)
            if(var$layout == "SINGLE") {
                r <- get_fastq_paths(fastqdir, var$run, cmd = "S")
                args <- c(refs, "--reads", r, "--threads", threads,
                          "--workdir", workdir, "--fastx",
                          "--aligned", paste0(workdir, "_rRNA"),
                          "--other", paste0(workdir, "_filt"))
                system2("sortmerna", args = args)
            } else if(var$layout == "PAIRED") {
                r1 <- get_fastq_paths(fastqdir, var$run, cmd = "P1")
                r2 <- get_fastq_paths(fastqdir, var$run, cmd = "P2")
                args <- c(refs, "--reads", r1, "--reads", r2,
                          "--threads", threads, 
                          "--workdir", workdir, "--fastx",
                          "--aligned", paste0(workdir, "_rRNA"),
                          "--other", paste0(workdir, "_filt"),
                          "--paired_in", "--out2")
                system2("sortmerna", args = args)
            } else {
                message("Layout information not available.")
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



