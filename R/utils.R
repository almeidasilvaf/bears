
#' Wrapper to extract variables for each Run or BioSample
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param index Row index.
#' 
#' @return A named list with:
#' \describe{
#'   \item{biosample}{BioSample accession.}
#'   \item{experiment}{Experiment accession.}
#'   \item{run}{Run accession.}
#'   \item{platform}{Sequencing platform.}
#'   \item{layout}{Library layout.}
#' }
#' @noRd
var2list <- function(sample_info, index=NULL) {
    biosample <- sample_info[index, "BioSample"]
    experiment <- sample_info[index, "Experiment"]
    run <- sample_info[index, "Run"]
    platform <- sample_info[index, "Instrument"]
    layout <- sample_info[index, "Layout"]
    
    res_list <- list(
        biosample = biosample, experiment = experiment,
        run = run, platform = platform, layout = layout
    )
    return(res_list)
}


#' Convert GFF file to BED
#'
#' @param gffpath Path to .gff file with genome annotation.
#' 
#' @return NULL
#' @export
#' @rdname gff2bed
#' @importFrom rtracklayer import export.bed
gff2bed <- function(gffpath=NULL) {
    gff <- rtracklayer::import(gffpath)
    gff$score <- as.numeric(rep(0, length(gff)))
    bedfile <- gsub(".gff", ".bed", gffpath)
    rtracklayer::export.bed(gff, bedfile)
    return(NULL)
}

#' Infer library strandedness
#'
#' @param mapping_passed Metadata of samples that passed mapping QC. This can be obtained with
#' \code{mapping_pass}.
#' @param bedpath Path to BED file. GFF files can be converted to BED with 
#' \code{gff2bed}.
#' @param mappingdir Directory where .bam files are stored.
#' 
#' @return A data frame with sample metadata as in mapping_passed, but with
#' an additional column named 'Orientation' containing library strandedness
#' for each BioProject.
#' @export
#' @rdname infer_strandedness
infer_strandedness <- function(mapping_passed = NULL,
                               bedpath = NULL,
                               mappingdir="results/04_read_mapping") {
    stranddir <- paste0(mappingdir, "/strandedness")
    # Create dir to store files with strandedness info
    if(!dir.exists(stranddir)) { dir.create(stranddir, recursive=TRUE) }
    selection <- mapping_passed[!duplicated(mapping_passed$BioProject), ]
    s <- lapply(seq_len(nrow(selection)), function(x) {
        bam <- paste0(mappingdir, "/", selection[x, "BioSample"],
                      "Aligned.sortedByCoord.out.bam")
        args <- c("-i", bam, "-r", bedpath, "-s 400000 >", 
                  paste0(stranddir, "/", selection[x, "BioSample"], ".txt"))
        system2("infer_experiment.py", args = args)
    })
    # Read output and get important info
    outfiles <- list.files(stranddir, full.names = TRUE)
    filelist <- lapply(outfiles, read.csv, header=FALSE, sep=" ", skip=3)
    strandedness <- unlist(lapply(filelist, function(x) {
        failed <- x$V7[1]
        explainedin <- x$V7[3]
        std <- sd(c(failed, explainedin))
        if(std < 0.1) {
            strand <- "unstranded"
        } else if(explainedin > failed) {
            strand <- "first"
        } else {
            strand <- "second"
        }
        return(strand)
    }))
    df <- data.frame(BioProject = selection$BioProject,
                     Orientation = strandedness)
    result <- merge(mapping_passed, df, by="BioProject", all.x=TRUE)
    return(result)
}
    