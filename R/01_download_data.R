
#' Wrapper to check if character vector is empty
#' 
#' @param vector Character vector to be inspected.
#' @param unique Logical indicating whether or not to pick only 
#' the first element if the character vector has 2 or more elements.
#' @return If the vector is empty, it is assigned NA. Otherwise, it remains
#' the same.
#' @noRd
check_empty <- function(vector, unique=TRUE) {
    if(unique) {
        vector <- vector[1]
    } else {
        vector <- vector
    }
    
    if(length(vector) == 0 | is.null(vector)) {
        vector <- NA
    }
    
    if(is.list(vector)) {
        vector <- NA
    }
    return(vector)
}


#' Create metadata data frame
#' 
#' @param res_list List containing variables to include in the metadata 
#' data frame.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{BioSample}{BioSample accession.}
#'   \item{Experiment}{Experiment accession (SRX*).}
#'   \item{Run}{Run accession (SRR*).}
#'   \item{Tissue}{Tissue from which RNA was extracted.}
#'   \item{Pubmed}{Pubmed ID of articles associated with the project, if any.}
#'   \item{BioProject}{Bioproject accession.}
#'   \item{Instrument}{Sequencing instrument}
#'   \item{Layout}{Library layout (single- or paired-end sequencing).}
#'   \item{Selection_method}{Library selection method.}
#'   \item{SRA_sample}{SRA sample accession (SRS*).}
#'   \item{SRA_study}{SRA study accession (SRP*).}
#'   \item{Treatment}{Treatment for the biosample.}
#'   \item{Cultivar}{Plant cultivar.}
#'   \item{Study_title}{Study title, if any.}
#'   \item{Study_abstract}{Study abstract, if any.}
#'   \item{Date}{Date of release.}
#'   \item{Origin}{Country of origin.}
#' }
#' @noRd
create_meta_df <- function(res_list) {
    nruns <- length(res_list$run)
    r <- res_list
    df <- data.frame(
        BioSample = rep(check_empty(r$biosample), nruns), 
        Experiment = rep(check_empty(r$experiment), nruns), 
        Run = check_empty(r$run, unique=FALSE),
        Tissue = rep(check_empty(r$tissue), nruns), 
        Pubmed = rep(check_empty(r$pubmed), nruns),
        BioProject = rep(check_empty(r$bioproject), nruns), 
        Instrument = rep(check_empty(r$instrument), nruns),
        Layout = rep(check_empty(r$layout), nruns), 
        Selection_method = rep(check_empty(r$selection), nruns),
        SRA_sample = rep(check_empty(r$srasample), nruns), 
        SRA_study = rep(check_empty(r$srastudy), nruns),
        Treatment = rep(check_empty(r$treatment), nruns), 
        Cultivar = rep(check_empty(r$cultivar), nruns),
        Study_title = rep(check_empty(r$title), nruns), 
        Study_abstract = rep(check_empty(r$abstract), nruns),
        Date = rep(as.character(check_empty(r$date))[1], nruns), 
        Origin = rep(check_empty(r$origin), nruns)
    )
    return(df)
}

#' Wrapper to create a data frame from efetch results in XML
#' 
#' @param id ID obtained from \code{rentrez::entrez_search()}.
#' 
#' @return A data frame with the following columns:
#' \describe{
#'   \item{BioSample}{BioSample accession.}
#'   \item{Experiment}{Experiment accession (SRX*).}
#'   \item{Run}{Run accession (SRR*).}
#'   \item{Tissue}{Tissue from which RNA was extracted.}
#'   \item{Pubmed}{Pubmed ID of articles associated with the project, if any.}
#'   \item{BioProject}{Bioproject accession.}
#'   \item{Instrument}{Sequencing instrument}
#'   \item{Layout}{Library layout (single- or paired-end sequencing).}
#'   \item{Selection_method}{Library selection method.}
#'   \item{SRA_sample}{SRA sample accession (SRS*).}
#'   \item{SRA_study}{SRA study accession (SRP*).}
#'   \item{Treatment}{Treatment for the biosample.}
#'   \item{Cultivar}{Plant cultivar.}
#'   \item{Study_title}{Study title, if any.}
#'   \item{Study_abstract}{Study abstract, if any.}
#'   \item{Date}{Date of release.}
#'   \item{Origin}{Country of origin.}
#' }
#' @noRd
#' @importFrom rentrez entrez_fetch
#' @importFrom XML xpathSApply xmlValue xmlChildren xmlAttrs
sra_xml2df <- function(id) {
    run_info <- rentrez::entrez_fetch(db="sra", id=id,
                                      rettype = "xml", parsed=TRUE)
    ext_id <- XML::xpathSApply(run_info, "//EXTERNAL_ID", XML::xmlValue)
    p_id <- XML::xpathSApply(run_info, "//IDENTIFIERS/PRIMARY_ID", XML::xmlValue)
    samp_at <- XML::xpathSApply(run_info, '//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE',
                                XML::xmlValue)
    biosample <- unique(ext_id[grep("SAM", ext_id)])
    bioproject <- unique(ext_id[grep("PRJ", ext_id)])
    experiment <- XML::xpathSApply(run_info, "//EXPERIMENT/IDENTIFIERS", 
                                   XML::xmlValue)
    run <- unique(p_id[grep(".+RR.*", p_id)])
    srasample <- unique(p_id[grep(".+RS.*", p_id)])
    srastudy <- unique(p_id[grep(".+RP.*", p_id)])
    tissue <- gsub("tissue", "", samp_at[grep("tissue", samp_at)])
    cultivar <- gsub("cultivar", "", samp_at[grep("cultivar", samp_at)])
    origin <- gsub("geo_loc_name", "", samp_at[grep("geo_loc_name", samp_at)])
    layout <- names(XML::xpathSApply(run_info, "//LIBRARY_LAYOUT", XML::xmlChildren))
    pubmed <- XML::xpathSApply(run_info, "//XREF_LINK", XML::xmlValue)
    pubmed <- gsub("pubmed", "", pubmed[grep("pubmed", pubmed)])
    selection <- XML::xpathSApply(run_info, "//LIBRARY_SELECTION", XML::xmlValue)
    instr <- XML::xpathSApply(run_info, "//PLATFORM", XML::xmlValue)
    title <- XML::xpathSApply(run_info, "//STUDY_TITLE", XML::xmlValue)
    abstract <- XML::xpathSApply(run_info, "//STUDY_ABSTRACT", XML::xmlValue)
    treatment <- gsub("treatment", "", samp_at[grep("treatment", samp_at)])
    date <- XML::xpathSApply(run_info, "//RUN", XML::xmlAttrs)
    if(methods::is(date, "matrix")) {
        date <- date["published", ]
    } else {
        date <- date[[1]]["published"]
    }
    
    res_list <- list(
        biosample = biosample, experiment = experiment, run = run, 
        tissue = tissue, pubmed = pubmed, bioproject = bioproject, 
        instrument = instr, layout = layout, selection = selection, 
        srasample = srasample, srastudy = srastudy, treatment = treatment,
        cultivar = cultivar, title = title, abstract = abstract, date = date,
        origin = origin
    )
    df <- create_meta_df(res_list)
    return(df)
} 

#' Search the SRA database and create a data frame of sample metadata
#' 
#' @param term Character with the search term, 
#' e.g. "Glycine max\[ORGN\] AND RNA-seq\[STRA\]".
#' @param retmax Numeric with the maximum number of hits returned 
#' by the search.
#' @param bp_param BiocParallel back-end to be used. 
#' Default: BiocParallel::SerialParam().
#' 
#' @return A data frame with the following columns:
#' \describe{
#'   \item{BioSample}{BioSample accession.}
#'   \item{Experiment}{Experiment accession (SRX*).}
#'   \item{Run}{Run accession (SRR*).}
#'   \item{Tissue}{Tissue from which RNA was extracted.}
#'   \item{Pubmed}{Pubmed ID of articles associated with the project, if any.}
#'   \item{BioProject}{Bioproject accession.}
#'   \item{Instrument}{Sequencing instrument}
#'   \item{Layout}{Library layout (single- or paired-end sequencing).}
#'   \item{Selection_method}{Library selection method.}
#'   \item{SRA_sample}{SRA sample accession (SRS*).}
#'   \item{SRA_study}{SRA study accession (SRP*).}
#'   \item{Treatment}{Treatment for the biosample.}
#'   \item{Cultivar}{Plant cultivar.}
#'   \item{Study_title}{Study title, if any.}
#'   \item{Study_abstract}{Study abstract, if any.}
#'   \item{Date}{Date of release.}
#'   \item{Origin}{Country of origin.}
#' }
#' @export
#' @rdname create_sample_info
#' @importFrom rentrez entrez_search
#' @importFrom BiocParallel bplapply SerialParam
#' @examples 
#' term <- "SAMN02422669[BSPL]"
#' df <- create_sample_info(term)
create_sample_info <- function(term, retmax=5000, 
                               bp_param = BiocParallel::SerialParam()) {
    search <- rentrez::entrez_search(db="sra", term=term,
                                     retmax = retmax, use_history = TRUE)
    final_list <- BiocParallel::bplapply(search$ids, sra_xml2df, 
                                         BPPARAM = bp_param)
    final_df <- Reduce(rbind, final_list)
    return(final_df)
}


#' Download FASTQ files
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' Default: results/01_FASTQ_files.
#' @param threads Number of threads to use. Default: 2.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
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
                           threads = 2,
                           envname = NULL, 
                           miniconda_path = NULL) {
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
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
fastq_exists <- function(sample_info = NULL, 
                         fastqdir = "results/01_FASTQ_files") {
    
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
    file_df$Run <- gsub("_.*|\\.fastq.*", "", file_df$Run)
    file_df <- file_df[!duplicated(file_df$Run), ]
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
        ftp_url <- utils::read.table(l, header=TRUE)$submitted_ftp
        Sys.sleep(1)
        ftp_url <- unlist(strsplit(ftp_url, ";"))
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
#' @examples 
#' data(sample_info)
#' fastqdir <- tempdir()
#' \donttest{
#' download_from_ena(sample_info, fastqdir)
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
        x <- download.file(urls[x], destfile = file, method = method)
    })
    
    df <- fastq_exists(sample_info, fastqdir)
    return(df)
     
}






