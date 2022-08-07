

#' Download FASTQ files
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' Default: results/01_FASTQ_files.
#' @param threads Number of threads to use. Default: 1.
#' 
#' @return A 2-column data frame with run accession in the first column
#' and status in the second column. If file is present, the status "OK"
#' is returned, otherwise NA is returned.
#' @export
#' @rdname download_fastq
#' @importFrom fs dir_delete
#' @examples
#' \donttest{
#' data(sample_info)
#' fastqdir <- tempdir()
#' if(sratoolkit_is_installed()) {
#'     download_fastq(sample_info, fastqdir)
#' }
#' }
download_fastq <- function(sample_info, 
                           fastqdir = "results/01_FASTQ_files",
                           threads = 1) {
    if(!sratoolkit_is_installed()) { stop("Unable to find SRAToolkit in PATH.") }

    if(!dir.exists(fastqdir)) { dir.create(fastqdir, recursive = TRUE) }
    d <- lapply(seq_len(nrow(sample_info)), function(x) {
        var <- var2list(sample_info, index = x)
        if(skip(var$platform)) {
            message("Skipping files...")
        } else {
            args <- c(var$run, "-e", threads, "-p --outdir", fastqdir)
            if(var$layout == "PAIRED") {
                args <- c(args, "--split-files")
            }
            system2("fasterq-dump", args = args)
        }
    })
    message("Compressing .fastq files...")
    system2("gzip", args = paste0(fastqdir, "/*.fastq"))
    flist <- list.files(path = fastqdir, pattern=".fastq.gz")
    flist <- unique(gsub("_[0-9].fastq.gz|.fastq.gz", "", flist))
    df_status <- data.frame(run = sample_info$Run)
    df_status$status <- ifelse(flist %in% df_status$run, "OK", NA)
    return(df_status)
}


#' Check if FASTQ files were properly downloaded
#'
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' Default: results/01_FASTQ_files.
#' @param collapse_pe Logical scalar indicating whether to collapse paired-end
#' files runs into a single run. If TRUE (default), 
#' files like SRR12345_1.fastq.gz and SRR12345_2.fastq.gz will be collapsed
#' to SRR12345.
#' 
#' @return Data frame with run accession in the first column, and
#' status in the second column. If the FASTQ file for a given run exists,
#' its status will be "OK", otherwise it will be NA.
#' @export
#' @rdname fastq_exists
#' @examples
#' data(sample_info)
#' fastqdir <- system.file("extdata", package = "bears")
#' fastq_exists(sample_info, fastqdir)
#' fastq_exists(sample_info, fastqdir, collapse_pe = FALSE)
fastq_exists <- function(sample_info = NULL, 
                         fastqdir = "results/01_FASTQ_files",
                         collapse_pe = TRUE) {
    
    files <- lapply(seq_len(nrow(sample_info)), function(x) {
        vars <- var2list(sample_info, index = x)
        if(vars$layout == "PAIRED") {
            file <- c(
                paste0(fastqdir, "/", vars$run, "_1.fastq.gz"),
                paste0(fastqdir, "/", vars$run, "_2.fastq.gz")
            )
        } else {
            file <- paste0(fastqdir, "/", vars$run, ".fastq.gz")
        }
        return(file)
    })
    file_df <- data.frame(Run = unlist(files), Status = NA)
    file_df$Status <- ifelse(file.exists(file_df$Run), "OK", NA)
    file_df$Run <- vapply(strsplit(file_df$Run, "/"), tail, n=1, character(1))
    file_df$Run <- gsub("\\.fastq.*", "", file_df$Run)
    if(collapse_pe) {
        file_df$Run <- gsub("_.*", "", file_df$Run)
        file_df <- file_df[!duplicated(file_df$Run), ]
    }
    return(file_df)
}


#' Get URL for each FASTQ file in the ENA's ftp repository via API
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' 
#' @return A character vector with the URL for each file.
#' @importFrom utils read.table
#' @noRd
get_url_ena_api <- function(sample_info = NULL) {
    
    base_url <- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession="
    runs <- sample_info$Run
    link <- unlist(lapply(runs, function(x) {
        l <- paste0(base_url, x, 
                    "&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&limit=0")
        
        ftp_url <- NULL
        try <- tryCatch({
            url <- utils::read.table(l, sep = "\t", header = TRUE)
            url <- url$fastq_ftp
            ftp_url <- unlist(strsplit(url, ";"))
            Sys.sleep(1)
        },
        error = function(e) {
            message("Could not find URL for run ", x)
            return(NULL)
        },
        warning = function(w) {
            message("Could not find URL for run ", x)
            return(NULL)
        })
        return(ftp_url)
    }))
    return(link)
}

#' Get URL for each FASTQ file in the ENA's FTP repository via iterations on possible links
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' 
#' @return A character vector with the URL for each file.
#' @noRd
get_url_ena_iterative <- function(sample_info = NULL) {
    
    base_url <- "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
    urls <- lapply(seq_len(nrow(sample_info)), function(x) {
        run <- sample_info$Run[x]
        layout <- sample_info$Layout[x]
        ext <- ".fastq.gz"
        subdir <- paste0(substr(run, 1, 6), "/")
        if(layout == "PAIRED") { ext <- c("_1.fastq.gz", "_2.fastq.gz") }
        file <- paste0(run, ext) 
        if(startsWith(run, "SRR")) {
            ssubdir <- paste0("00", substr(run, nchar(run), nchar(run)), "/")
            url <- paste0(base_url, subdir, ssubdir, paste0(run, "/"), file)
        } else {
            url <- paste0(base_url, subdir, paste0(run, "/"), file)
        }
        check <- vapply(url, valid_url, logical(1))
        if(any(check == FALSE)) {
            ssubdir <- paste0("0", substr(run, nchar(run)-1, nchar(run)), "/")
            url <- paste0(base_url, subdir, ssubdir, paste0(run, "/"), file)
            check <- vapply(url, valid_url, logical(1))
            
            if(any(check == FALSE)) {
                url <- gsub("/[0-9][0-9][0-9]/", "/", url)
                check <- vapply(url, valid_url, logical(1))
                if(any(check == FALSE)) {
                    message("Could not find URL for run ", run)
                    url <- NULL
                }
            }
        }
        return(url)
    })
    urls <- urls[!vapply(urls, is.null, logical(1))]
    urls <- unlist(urls)
    return(urls)
}


#' Get URL for each file in the ENA's FTP repository
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param link_from Method to extract the URL to each FASTQ file in the ENA's
#' ftp repository. One of 'api' or 'iterative'. Default: 'api'.
#' 
#' @return A character vector with the URL for each accession.
#' @export
#' @rdname get_url_ena
#' @examples 
#' data(sample_info)
#' get_url_ena(sample_info)
get_url_ena <- function(sample_info = NULL, link_from = "api") {
    
    check_internet <- valid_url("https://google.com")
    if(!check_internet) { stop("You have an internet connection problem.") }
    
    if(link_from == "api") {
        urls <- get_url_ena_api(sample_info)
    } else if(link_from == "iterative") {
        urls <- get_url_ena_iterative(sample_info)
    } else {
        stop("Invalid 'link_from'. Choose one of 'api' or 'iterative'.")
    }
    return(urls)
}

#' Download FASTQ files from ENA's FTP
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param urls Character vector returned by \code{get_url_ena()} with
#' the URLs to each file in the ENA's FTP repository. If NULL, this function
#' will run \code{get_url_ena()} to get the URLs before downloading.
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' Default: results/01_FASTQ_files.
#' @param method Method to be used for downloading files. One of "internal",
#' "libcurl", "wget", "libcurl", "curl", "wininet" (Windows only), or "auto".
#' @param link_from Method to extract the URL to each FASTQ file in the ENA's
#' ftp repository. One of 'api' or 'iterative'. Default: 'api'.
#' 
#' @return A data frame as returned by \code{fastq_exists}.
#' @rdname download_from_ena
#' @export
#' @importFrom utils download.file tail
#' @importFrom downloader download
#' @examples 
#' data(sample_info)
#' fastqdir <- tempdir()
#' \donttest{
#' download_from_ena(sample_info, fastqdir = fastqdir)
#' }
download_from_ena <- function(sample_info = NULL, 
                              urls = NULL,
                              fastqdir = "results/01_FASTQ_files", 
                              method = "auto",
                              link_from = "api") {
    if(missing(method)) 
        method <- ifelse(!is.null(getOption("download.file.method")), 
                         getOption("download.file.method"), "auto")
    
    if(is.null(urls)) {
        urls <- get_url_ena(sample_info, link_from = link_from)
    }
    d <- lapply(seq_along(urls), function(x) {
        message("Downloading file ", urls[x])
        file <- vapply(strsplit(urls[x], "/"), tail, n=1, character(1))
        file <- paste0(fastqdir, "/", file)
        # Try to download file: if it doesn't work, delete intermediate file
        x <- tryCatch({
            downloader::download(urls[x], destfile = file, method = method)
            res <- TRUE
        },
        error = function(e) {
            message("Could not download file ", urls[x])
            return(FALSE)
        },
        warning = function(w) {
            message("Could not download file ", urls[x])
            return(FALSE)
        })
        if(!x) { unlink(file) }
    })
    
    df <- fastq_exists(sample_info, fastqdir)
    return(df)
     
}


#' Check integrity of downloaded FASTQ files
#' 
#' This file detects FASTQ files that have not been properly downloaded
#' and flags them, indicating the problem in a "Issue" column. It is an
#' alternative to md5sum to check file integrity.
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info()}.
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' Default: results/01_FASTQ_files. 
#' @param read_count A 2-column data frame with the number of reads for
#' each run as reported in SRA, which can be obtained with the
#' function \code{get_read_count()}.
#' @param verbose Logical scalar indicating whether or not to print 
#' log messages. Default: FALSE.
#'
#' @return A 3-column data frame with the following variables:
#' \describe{
#'   \item{Run}{Character vector of run acessions.}
#'   \item{Status}{Character vector of file integrity status. 
#'   If there is nothing wrong with the file, it will be "OK", 
#'   and NA otherwise.}
#'   \item{Issue}{Numeric vector with code to issue in each file, if there is
#'   any. See issue legend below.}
#' }
#' @details
#' This function looks for and flags problematic downloads.
#' Run accessions that have been properly downloaded will have "OK" in 
#' the \strong{Status} variable. If there is any issue with the run, 
#' the \strong{Status} variable will be NA, with details on the problem
#' in the \strong{Issue} column. Issues are labeled numbers from 1 to 3,
#' which mean:
#' * 1: File was not downloaded.
#' * 2: For paired-end reads, only forward or reverse file was downloaded.
#' * 3: Less reads than the expected. The expected number of reads for each
#' file is available in SRA, and it can be obtained 
#' with \code{get_read_count()}.
#' 
#' @export
#' @importFrom ShortRead countFastq
#' @rdname check_downloads
#' @examples
#' data(sample_info)
#' fastqdir <- system.file("extdata", package = "bears") 
#' # Create a fake data frame of read count
#' read_count <- data.frame(Run = "SRR1039508", Reads = 7097)
#' # Check download
#' check_downloads(sample_info, fastqdir, read_count)
check_downloads <- function(sample_info = NULL, 
                            fastqdir = "results/01_FASTQ_files",
                            read_count,
                            verbose = FALSE) {
    
    # Get data frame of downloaded files
    downloaded <- fastq_exists(sample_info, fastqdir, collapse_pe = FALSE)
    downloaded$CRun <- gsub("_[0-9]", "", downloaded$Run)
    
    # Create list of runs
    run_list <- split(downloaded, downloaded$CRun)
    
    # Define wrapper function to check if number of reads match the expected
    check_nreads <- function(fastqdir, run, read_count) {
        file <- list.files(fastqdir, pattern = run[["CRun"]], full.names = TRUE)
        file <- file[endsWith(file, ".fastq.gz")]
        run$nreads <- ShortRead::countFastq(file)$records
        run$nref <- read_count$Reads[read_count$Run %in% run$CRun]
        return(identical(run$nreads, run$nref))
    }
    
    # Add status column
    run_status <- Reduce(rbind, lapply(run_list, function(x) {
        if(verbose) {
            message("Working on run ", unique(x$CRun))
        }
        x$Issue <- NA
        
        ## Single-end file
        if(nrow(x) == 1) { 
            if(is.na(x$Status)) {
                x$Issue <- 1 
            } else {
                ### Check if number of reads match the reported in SRA
                check_readcount <- check_nreads(fastqdir, x, read_count) 
                if(!check_readcount) {
                    x$Status <- NA
                    x$Issue <- 3
                }
            }
            
            ## Paired-end file
        } else {
            ### Check if any file is missing. None should be.
            na_count <- sum(is.na(x$Status))
            
            if(na_count == 2) { # both forward and reverse not present
                x$Status <- NA
                x$Issue <- 1
            } else if(na_count == 1) { # either forward or reverse not present
                x$Status <- NA
                x$Issue <- 2
            } else { # Both files were downloaded
                ### Check if number of reads match the reported in SRA
                check_readcount <- check_nreads(fastqdir, x, read_count) 
                if(!check_readcount) {
                    x$Status <- NA
                    x$Issue <- 3
                }
            }
        }
        final_status <- x[1, c("CRun", "Status", "Issue")]
        names(final_status) <- c("Run", "Status", "Issue")
        if(verbose & !is.na(final_status$Issue)) {
            message("Issue ", final_status$Issue, " in run ", final_status$Run)
        }
        return(final_status)
    }))
    return(run_status)
}
