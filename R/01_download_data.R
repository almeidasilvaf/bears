

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
        
        failed <- function(run) {
            message("Could not find URL for run ", run)
            return(NULL)
        }
        
        ftp_url <- NULL
        try <- tryCatch({
            url <- utils::read.table(l, sep = "\t", header = TRUE)
            url <- url$fastq_ftp
            ftp_url <- unlist(strsplit(url, ";"))
            Sys.sleep(1)
        },
        error = function(e) { failed(x) },
        warning = function(w) { failed(x) }
        )
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
    
    # Handle paired-end reads with 3 urls - keep only _1 and _2.fastq.gz
    todelete <- unlist(lapply(sample_info$Run, function(x) {
        count <- sum(grepl(x, urls))
        delete <- "nothing"
        if(count > 2) {
            delete <- x
        }
        return(delete)
    }))
    urls <- urls[!grepl(paste0(todelete, ".fastq.gz"), urls)]

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

#' Check file integrity with md5sum
#'
#' @param run_accessions Character vector of run accessions.
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' Default: results/01_FASTQ_files. 
#' 
#' @return A data frame with variables \strong{Run} and \strong{Status} with
#' run accession and integrity status, respectively. 
#' 
#' @export
#' @rdname check_md5
#' @examples
#' urls <- c(
#'     "ftp.sra.ebi.ac.uk/vol1/fastq/SRR926/SRR926397/SRR926397_1.fastq.gz",
#'     "ftp.sra.ebi.ac.uk/vol1/fastq/SRR926/SRR926397/SRR926397_2.fastq.gz"
#' )
#' sample_info <- data.frame(
#'     BioSample = "SAMN01924555",
#'     Experiment = "SRX245306",
#'     Run = "SRR926397",
#'     BioProject = "PRJNA190191", Instrument = "Illumina HiSeq 2000", 
#'     Layout = "PAIRED"
#' )
#' fastqdir <- tempdir()
#' d <- download_from_ena(
#'     sample_info, urls = urls, fastqdir = fastqdir, method = "libcurl"
#' )
#' 
#' # Check MD5
#' run_accessions <- sample_info$Run
#' check_md5(run_accessions, fastqdir)
check_md5 <- function(run_accessions = NULL, 
                      fastqdir = "results/01_FASTQ_files") {
    
    base_url <- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession="
    check <- Reduce(rbind, lapply(run_accessions, function(x) {
        
        l <- paste0(base_url, x, "&result=read_run&format=tsv&limit=0")
        
        md5_status <- data.frame(Run = x, Status = NA)
        try <- tryCatch({
            info <- utils::read.table(l, sep = "\t", header = TRUE)
            md5 <- data.frame(
                Run = unlist(strsplit(info$fastq_ftp, ";")),
                md5 = unlist(strsplit(info$fastq_md5, ";"))
            )
            # Remove cases with 3 FASTQ files for paired-end reads
            paired_3 <- sum(
                grepl(paste0(x, ".fastq.gz"), md5$Run),
                grepl(paste0(x, "_1.fastq.gz"), md5$Run)
            )
            if(paired_3 != 1) {
                md5 <- md5[!grepl(paste0(x, ".fastq.gz"), md5$Run), ]
            }
            
            md5$md5_downloaded <- tools::md5sum(
                file.path(fastqdir, basename(md5$Run))
            )
            md5$Status <- identical(md5$md5, md5$md5_downloaded)
            
            # Get status
            md5_status$Status <- ifelse(all(md5$Status == TRUE), TRUE, FALSE)
        },
        error = function(e) {
            message("Could not find URL for run ", x)
            return(NULL)
        },
        warning = function(w) {
            message("Could not find URL for run ", x)
            return(NULL)
        })
        
        return(md5_status)
    }))
    
    return(check)
}

