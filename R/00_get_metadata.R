
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
    if(length(res_list$run) != 0) {
        nruns <- length(res_list$run)
    } else {
        nruns <- length(res_list$experiment)
    }
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
    rawdate <- XML::xpathSApply(run_info, "//RUN", XML::xmlAttrs)
    date <- NA
    if(methods::is(rawdate, "matrix")) {
        date <- rawdate["published", ]
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
#' term <- "SAMN02422669[BSPL]"
#' df <- create_sample_info(term)
create_sample_info <- function(term, retmax=5000) {
    search <- rentrez::entrez_search(db="sra", term=term,
                                     retmax = retmax, use_history = TRUE)
    final_list <- lapply(search$ids, sra_xml2df)
    final_df <- Reduce(rbind, final_list)
    return(final_df)
}

