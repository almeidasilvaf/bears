
#' Wrapper to skip reads
#' 
#' Files are skipped if they do not exist where they should or if they are
#' SOLiD/PacBio reads.
#' 
#' @param platform Sequencing platform.
#' @param path Path to file to be tested.
#' 
#' @return Logical indicating whether to skip reads or not.
#' @noRd
skip <- function(platform = NULL, path = NULL) {
    p_skip <- FALSE
    if(grepl("SOLiD|PacBio", platform)) {
        p_skip <- TRUE
    }
    if(is.null(path)) {
        f_skip <- FALSE
    } else {
        f_skip <- !file.exists(path)
    }
    final <- FALSE
    if(p_skip | f_skip) {
        final <- TRUE
    }
    return(final)
}

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
#' @return A 2-column data frame with path to BED file in the first column
#' and status in the second column, with "OK" if the BED file was successfully
#' created from the GFF/GTF file, and NA otherwise.
#' @export
#' @rdname gff2bed
#' @importFrom rtracklayer import export.bed
#' @examples 
#' gff_path <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf", 
#'                          package="bears")
#' gffdir <- tempdir()
#' file.copy(from = gff_path, to=gffdir)
#' gff_file <- list.files(gffdir, full.names=TRUE, pattern=".gtf")
#' gff2bed(gff_file)
gff2bed <- function(gffpath=NULL) {
    gff <- rtracklayer::import(gffpath)
    gff$score <- as.numeric(rep(0, length(gff)))
    bedfile <- gsub(".gff|.gtf|.gff3", ".bed", gffpath)
    rtracklayer::export.bed(gff, bedfile)
    status <- "OK"
    if(!file.exists(bedfile)) { status <- NA }
    df_status <- data.frame(bed_path = bedfile,
                            status = status)
    return(df_status)
}

#' Infer library strandedness
#'
#' @param mapping_passed Metadata of samples that passed mapping QC. 
#' This can be obtained with \code{mapping_pass}.
#' @param bedpath Path to BED file. GFF files can be converted to BED with 
#' \code{gff2bed}.
#' @param mappingdir Directory where .bam files are stored.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#' 
#' @return A data frame with sample metadata as in mapping_passed, but with
#' an additional column named 'Orientation' containing library strandedness
#' for each BioProject.
#' @export
#' @rdname infer_strandedness
#' @importFrom utils read.csv
#' @importFrom stats sd
#' @examples
#' data(sample_info)
#' mapping_passed <- sample_info
#' bedpath <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.bed", 
#'                         package="bears")
#' mappingdir <- system.file("extdata", package="bears")
#' if(rseqc_is_installed()) {
#'     strandedness <- infer_strandedness(mapping_passed, bedpath, mappingdir)
#' }
#' 
infer_strandedness <- function(mapping_passed = NULL,
                               bedpath = NULL,
                               mappingdir="results/04_read_mapping",
                               envname = NULL,
                               miniconda_path = NULL) {
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!rseqc_is_installed()) { stop("Unable to find RSeQC in PATH.") }
    
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
        std <- stats::sd(c(failed, explainedin))
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


#' Translate library orientation terminology for each program
#' 
#' @param orientation Library orientation, available in 
#' the column "Orientation" from the output of \code{infer_strandedness}.
#' @param layout Library layout, available in the column "Layout" from the
#' output of \code{create_sample_info}.
#' 
#' @return A list with the following elements:
#' \describe{
#'   \item{salmon}{salmon library information.}
#'   \item{fcounts}{featureCounts library information.}
#'   \item{kallisto}{kallisto library information.}
#'   \item{stringtie}{StringTie library information.}
#' }
#' @rdname translate_strandedness
#' @export
#' @examples 
#' data(sample_info)
#' orientation <- sample_info$Orientation
#' layout <- sample_info$Layout
#' strandedness <- translate_strandedness(orientation, layout)
translate_strandedness <- function(orientation = NULL, layout = NULL) {
    if(orientation == "unstranded") {
        salmon <- "U"
        fcounts <- 0
        kallisto <- ""
        stringtie <- ""
    } else if(orientation == "first") {
        salmon <- "SR"
        fcounts <- 2
        kallisto <- "--rf-stranded"
        stringtie <- "--rf"
    } else if(orientation == "second") {
        salmon <- "SF"
        fcounts <- 1
        kallisto <- "--fr-stranded"
        stringtie <- "--fr"
    } else {
        salmon <- "unknown"
        fcounts <- "unknown"
        kallisto <- "unknown"
        stringtie <- "unknown"
    }
    if(layout == "PAIRED") {
        salmon <- paste0("I", salmon)
    }
    strandedness <- list(salmon = salmon, 
                         fcounts = fcounts, 
                         kallisto = kallisto,
                         stringtie = stringtie)
    return(strandedness)
}




 