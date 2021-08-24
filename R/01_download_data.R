
#' Wrapper to check if character vector is empty
#' 
#' @param vector Character vector to be inspected.
#' 
#' @return If the vector is empty, it is assigned NA. Otherwise, it remains
#' the same.
#' @noRd
check_empty <- function(vector) {
    if(length(vector) == 0 | is.null(vector)) {
        vector <- NA
    }
    return(vector)
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
    date <- XML::xpathSApply(run_info, "//RUN", XML::xmlAttrs)["published",]
    
    df <- data.frame(
        BioSample = check_empty(biosample), 
        Experiment = check_empty(experiment), Run = check_empty(run),
        Tissue = check_empty(tissue), Pubmed = check_empty(pubmed),
        BioProject = check_empty(bioproject), Instrument = check_empty(instr),
        Layout = check_empty(layout), Selection_method = check_empty(selection),
        SRA_sample = check_empty(srasample), SRA_study = check_empty(srastudy),
        Treatment = check_empty(treatment), Cultivar = check_empty(cultivar),
        Study_title = check_empty(title), Study_abstract = check_empty(abstract),
        Date = as.character(check_empty(date)), Origin = check_empty(origin)
    )
    return(df)
} 

#' Search the SRA database and create a data frame of sample metadata
#' 
#' @param term Character with the search term, 
#' e.g. "Glycine max\[ORGN\] AND RNA-seq\[STRA\]".
#' @param retmax Numeric with the maximum number of hits returned 
#' by the search.
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
#' @examples 
#' term <- "PRJNA80173[GPRJ]"
#' df <- create_sample_info(term)
create_sample_info <- function(term, retmax=5000) {
    search <- rentrez::entrez_search(db="sra", term=term,
                                     retmax = retmax, use_history = TRUE)
    final_df <- Reduce(rbind, lapply(search$ids, sra_xml2df))
    return(final_df)
}


#' Download FASTQ files
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param sradir Path to the directory where .sra files will be temporarily
#' stored. Default: results/00_SRA_files.
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' Default: results/01_FASTQ_files.
#' @param soliddir Path to the directory where .fastq files for SOLiD data
#' will be stored. Default: results/01_SOLiD_dir.
#' 
#' @return NULL, with .sra files in the directory specified in sradir, and 
#' fastq files in the directories specified in fastqdir and soliddir.
#' @export
#' @rdname download_fastq
#' @importFrom fs dir_delete
#' @examples
#' \donttest{
#' data(sample_info)
#' download_fastq(sample_info, threads = 20)
#' }
download_fastq <- function(sample_info, 
                           fastqdir="results/01_FASTQ_files",
                           sradir="results/00_SRA_files", 
                           soliddir="results/01_SOLiD_dir", 
                           threads = 6) {
    if(!dir.exists(fastqdir)) { dir.create(fastqdir, recursive = TRUE) }
    d <- lapply(seq_len(nrow(sample_info)), function(x) {
        var <- var2list(sample_info, index = x)
        if(grepl("SOLiD", var$platform)) {
            if(!dir.exists(soliddir)) { dir.create(soliddir, recursive = TRUE) }
            if(!dir.exists(sradir)) { dir.create(sradir, recursive = TRUE) }
            args <- c("progress 3 -O", sradir, var$run)
            system2("prefetch", args = args)
            
            system2("abi-dump", args = c(
                "--outdir", soliddir, paste0(sradir, "/", var$run, ".sra")  
            ))
            fs::dir_delete(sradir)
        } else {
            args <- c(var$run, "-e", threads, "-p --outdir", fastqdir)
            if(var$layout == "PAIRED") {
                args <- c(args, "--split-files")
            }
            system2("fasterq-dump", args = args)
        }
    })
    message("Compressing .fastq files...")
    system2("gzip", args = paste0(fastqdir, "/*"))
    return(NULL)
}


