
#' Index genome for STAR alignment
#'
#' @param genome_path Path to the genome .fasta file.
#' @param gff_path Path to the .gff/.gtf file with annotations.
#' @param mappingdir Path to the directory where read mapping files (.bam) will
#' be stored.
#' @param indexdir Directory where the STAR genome index files will be stored.
#' Default: results/04_read_mapping/genomeIndex.
#' @param threads Number of threads for STAR aligner. Default: 2.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#' 
#' @return A 2-column data frame with path to index in the first column and
#' index build status in the second column, with "OK" if the index was 
#' built successfully and NA otherwise.
#' @export
#' @rdname star_genome_index
#' @importFrom tools file_ext
#' @examples 
#' \donttest{
#' genome_path <- system.file("extdata", "Hsapiens_GRCh37.75_subset.fa", 
#'                             package="bears")
#' gff_path <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf", 
#'                          package="bears")
#' mappingdir <- tempdir()
#' indexdir <- tempdir()
#' threads <- 2
#' if(star_is_installed()) {
#'     star_genome_index(genome_path, gff_path, mapping_dir, indexdir)
#' }
#' }
star_genome_index <- function(genome_path = NULL,
                              gff_path = NULL,
                              mappingdir = "results/04_read_mapping",
                              indexdir = "results/04_read_mapping/genomeIndex",
                              threads = 2,
                              envname = NULL,
                              miniconda_path = NULL) {
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!star_is_installed()) { stop("Unable to find STAR in PATH.") }

    if(!dir.exists(indexdir)) { dir.create(indexdir, recursive = TRUE) }
    exonparent <- "Parent"
    if(tools::file_ext(gff_path) == "gtf") { exonparent <- "transcript_id" }
    args <- c("--runThreadN", threads, "--runMode genomeGenerate",
              "--genomeDir", indexdir, "--genomeFastaFiles", genome_path,
              "--sjdbGTFfile", gff_path, 
              "--sjdbGTFtagExonParentTranscript", exonparent, 
              "--sjdbOverhang 25")
    system2("STAR", args = args)
    status <- "OK"
    if(length(dir(indexdir)) == 0) { status <- NA }
    df_status <- data.frame(index_path = indexdir,
                            status = status)
    return(df_status)
}


#' Wrapper to create a list of paths to reads to pass to STAR
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param filtdir Path to the directory where filtered reads will be stored.
#' Default: results/03_filtered_FASTQ.
#' 
#' @return A list with paths to reads in a suitable format for STAR input.
#' @noRd
star_reads <- function(sample_info, 
                       filtdir = "results/03_filtered_FASTQ") {
    single <- sample_info[sample_info$Layout == "SINGLE", ]
    paired <- sample_info[sample_info$Layout == "PAIRED", ]
    slist <- split(single, single$BioSample)
    plist <- split(paired, paired$BioSample)
    
    slist2 <- lapply(slist, function(x) {
        x$Run <- paste0(filtdir, "/", x$Run, ".fastq.gz")
        if(nrow(x) > 1) {
            y <- paste0(x$Run, collapse = ",")
        } else {
            y <- x$Run
        }
        return(y)
    })
    
    plist2 <- lapply(plist, function(x) {
        x$Run <- paste0(filtdir, "/", x$Run)
        if(nrow(x) == 1) {
            p1 <- paste0(x$Run, "_1.fastq.gz")
            p2 <- paste0(x$Run, "_2.fastq.gz")
            y <- paste(p1, p2, sep = " ")
        } else {
            runs <- x$Run
            p1 <- vapply(runs, function(x) paste0(x, "_1.fastq.gz"), character(1))
            p1 <- paste(p1, collapse = ",")
            p2 <- vapply(runs, function(x) paste0(x, "_2.fastq.gz"), character(1))
            p2 <- paste(p2, collapse = ",")
            y <- paste(p1, p2, sep = " ")
        }
        return(y)
    })
    final_list <- c(slist2, plist2)
    return(final_list)
}


#' Align reads to a reference genome using STAR
#'
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param filtdir Path to the directory where filtered reads will be stored.
#' Default: results/03_filtered_FASTQ.
#' @param fastqc_table Data frame of summary statistics for FastQC as returned
#' by \code{multiqc()}.
#' @param mappingdir Path to the directory where read mapping files (.bam) will
#' be stored.
#' @param indexdir Directory where the STAR genome index files will be stored.
#' Default: results/04_read_mapping/genomeIndex.
#' @param gff_path Path to the .gff/.gtf file with annotations.
#' @param threads Number of threads for STAR aligner. Default: 2.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#'
#' @return A 2-column data frame with BioSample IDs in the first column
#' and STAR running status in the second column, with "OK" if reads 
#' were mapped (.bam files were created) and NA if STAR failed to map reads.
#' @importFrom tools file_ext
#' @export
#' @rdname star_align
#' @examples
#' \donttest{
#' data(sample_info)
#' data(fastqc_table)
#' genome_path <- system.file("extdata", "Hsapiens_GRCh37.75_subset.fa", 
#'                             package="bears")
#' gff_path <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf", 
#'                          package="bears")
#' mappingdir <- tempdir()
#' indexdir <- tempdir()
#' threads <- 2
#' filtdir <- system.file("extdata", package="bears")
#' if(star_is_installed()) {
#'     star_genome_index(genome_path, gff_path, mapping_dir, indexdir)
#'     star_align(sample_info, filtdir, fastqc_table, mappingdir, 
#'                indexdir, gff_path, threads)
#' }
#' }
star_align <- function(sample_info = NULL,
                       filtdir = "results/03_filtered_FASTQ",
                       fastqc_table = NULL,
                       mappingdir = "results/04_read_mapping",
                       indexdir = "results/04_read_mapping/genomeIndex",
                       gff_path = NULL,
                       threads = 2,
                       envname = NULL,
                       miniconda_path = NULL) {
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!star_is_installed()) { stop("Unable to find STAR in PATH.") }
    
    exonparent <- "Parent"
    if(tools::file_ext(gff_path) == "gtf") { exonparent <- "transcript_id" }
    fastqc_table$Sample <- gsub("_.*", "", fastqc_table$Sample)
    qc_table <- fastqc_table[!duplicated(fastqc_table$Sample), 
                             c("Sample", "Sequence.length")]
    sample_info2 <- merge(sample_info[!duplicated(sample_info$BioSample),],
                          qc_table, by.x="Run", "Sample")
    reads <- star_reads(sample_info, filtdir)
    aln <- lapply(seq_len(nrow(sample_info2)), function(x) {
        platform <- sample_info2[x, "Instrument"]
        if(grepl("SOLiD|PacBio", platform)) {
            message("Skipping PacBio/SOLiD reads...")
        } else {
            prefix <- paste0(mappingdir, "/", sample_info2[x, "BioSample"])
            args <- c("--runThreadN", threads, "--genomeDir", indexdir, 
                      "--readFilesIn", reads[[x]], "--readFilesCommand zcat",
                      "--outFileNamePrefix", prefix, 
                      "--outSAMtype BAM SortedByCoordinate",
                      "--outSAMstrandField intronMotif", 
                      "--sjdbGTFfile", gff_path, 
                      "--sjdbGTFtagExonParentTranscript", exonparent,
                      "--sjdbOverhang 25")
            system2("STAR", args = args)
        }
    })
    flist <- list.files(path = mappingdir, pattern = ".bam")
    flist <- gsub("Aligned.*", "", flist)
    df_status <- data.frame(sample = sample_info$BioSample)
    df_status$status <- ifelse(flist %in% df_status$sample, "OK", NA)
    return(df_status)
}


#' Keep only samples that passed minimum requirements in STAR alignment
#' 
#' @param mapping_qc Data frame with read mapping statistics, as generated by
#' \code{multiqc()}.
#' @param sample_info Data frame with sample metadata as generated by 
#' \code{create_sample_info}.
#' 
#' @details Samples are excluded if i. more than 50 percent of the reads
#' fail to map, or ii. more than 40 percent of the reads fail to uniquely map.
#' 
#' @return Data frame with metadata of samples that passed alignment QC.
#' @export
#' @rdname mapping_pass
#' @examples
#' data(sample_info)
#' data(mapping_qc)
#' mapping_passed <- mapping_pass(mapping_qc, sample_info)
mapping_pass <- function(mapping_qc = NULL,
                         sample_info = NULL) {
    filt <- mapping_qc
    filt$uniqueperc <- filt$uniquely_mapped_percent
    filt$unmapperc <- filt$unmapped_tooshort_percent + 
        filt$unmapped_other_percent + 
        filt$unmapped_mismatches_percent
    
    filt <- filt[filt$uniqueperc >= 50 & filt$unmapperc < 40, ]
    sinfo <- sample_info[sample_info$BioSample %in% filt$Sample, ]
    return(sinfo)
}






