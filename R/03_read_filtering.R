
#' Run Trimmomatic on reads that failed on FastQC check
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param fastqc_table Data frame of summary statistics for FastQC as returned
#' by \code{multiqc()}.
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' Default: results/01_FASTQ_files.
#' @param filtdir Path to the directory where filtered .fastq files will
#' be stored. Default: results/03_filtered_FASTQ.
#' 
#' @export
#' @rdname trim_reads
#' @importFrom fs file_move file_delete
#' @importFrom utils download.file
#' @return NULL
trim_reads <- function(sample_info = NULL,
                       fastqc_table = NULL,
                       fastqdir = "results/01_FASTQ_files",
                       filtdir = "results/03_filtered_FASTQ",
                       adapterdir = "results/03_filtered_FASTQ/adapters") {
    if(!dir.exists(filtdir)) { dir.create(filtdir, recursive = TRUE) }
    if(!dir.exists(adapterdir)) {
        dir.create(adapterdir, recursive = TRUE)
        pe_ad <- paste0(adapterdir, "/PE_adap.fa")
        se_ad <- paste0(adapterdir, "/SE_adap.fa")
        download.file("https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa",
                      destfile = pe_ad)
        download.file("https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-SE.fa",
                      destfile = se_ad)
    }
    failed <- fastqc_table[fastqc_table$per_base_sequence_quality == "fail", 1]
    failed <- unique(gsub("_1$|_2$", "", failed))
    sample_info <- sample_info[sample_info$Run %in% failed, ]
    if(nrow(sample_info) > 0) {
        t <- lapply(seq_len(nrow(sample_info)), function(x) {
            var <- var2list(sample_info, index = x)
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
    nonfilt <- 
    return(NULL)
}


#' Wrapper to download standard rRNA databases from SortMeRNA
#' 
#' @param filtdir Path to the directory where filtered reads will be stored.
#' Default: results/03_filtered_FASTQ
#' 
#' @return "--ref" parameters for SortMeRNA.
#' @noRd
download_rrna <- function(filtdir = "results/03_filtered_FASTQ") {
    outdir <- paste0(filtdir, "/rRNAdbs/")
    # Create paths to downloaded files
    fivep8s <- paste0(outdir, "rfam_5.8s.fasta")
    fives <- paste0(outdir, "rfam_5s.fasta")
    arc_16s <- paste0(outdir, "arc_16s.fasta")
    arc_23s <- paste0(outdir, "arc_23s.fasta")
    bac_16s <- paste0(outdir, "bac_16s.fasta")
    bac_23s <- paste0(outdir, "bac_23s.fasta")
    euk_18s <- paste0(outdir, "euk_18s.fasta")
    euk_28s <- paste0(outdir, "euk_28s.fasta")
    
    # Download files
    files <- c(fivep8s, fives, arc_16s, arc_23s, 
               bac_16s, bac_23s, euk_18s, euk_28s)
    if(sum(file.exists(files)) != 8) {
        download.file("https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/rfam-5.8s-database-id98.fasta",
                      destfile = fivep8s)
        download.file("https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/rfam-5s-database-id98.fasta",
                      destfile = fives)
        download.file("https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-arc-16s-id95.fasta",
                      destfile = arc_16s)
        download.file("https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-arc-23s-id98.fasta",
                      destfile = arc_23s)
        download.file("https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-bac-16s-id90.fasta",
                      destfile = bac_16s)
        download.file("https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-bac-23s-id98.fasta",
                      destfile = bac_23s)
        download.file("https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-euk-18s-id95.fasta",
                      destfile = euk_18s)
        download.file("https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-euk-28s-id98.fasta",
                      destfile = euk_28s)
    }
    refs <- c("--ref", fivep8s, "--ref", fives, 
              "--ref", arc_16s, "--ref", arc_23s,
              "--ref", bac_16s, "--ref", bac_23s,
              "--ref", euk_18s, "--ref", euk_28s)
    return(refs)
}


#' Wrapper to clean directory after running SortMeRna
#' 
#' @param filtdir Path to the directory where filtered reads will be stored.
#' Default: results/03_filtered_FASTQ.
#' 
#' @return NULL
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


#' Remove rRNA sequences from .fastq files with SortMeRna
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param fastqdir Path to the directory where .fastq files are stored.
#' @param filtdir Path to the directory where filtered reads will be stored.
#' Default: results/03_filtered_FASTQ.
#' @param threads Number of threads for SortMeRna. Default: 8.
#' 
#' @export
#' @rdname remove_rrna
#' @importFrom fs dir_delete
remove_rrna <- function(sample_info,
                        fastqdir = "results/01_FASTQ_files",
                        filtdir = "results/03_filtered_FASTQ",
                        threads = 8) {
    if(!dir.exists(filtdir)) { dir.create(filtdir, recursive = TRUE) }
    dbdir <- paste0(filtdir, "/rRNAdbs")
    if(!dir.exists(dbdir)) { dir.create(dbdir, recursive = TRUE) }
    refs <- download_rrna(filtdir)
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
    return(NULL)
}



