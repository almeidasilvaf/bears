

#' Index genome for STAR alignment
#'
#' @param genome_path Path to the genome .fasta file.
#' @param gff_path Path to the .gff/.gtf file with annotations.
#' @param mappingdir Path to the directory where read mapping files (.bam) will
#' be stored.
#' @param threads Number of threads for STAR aligner.
#' 
#' @return NULL
#' @export
#' @rdname star_genome_index
#' @importFrom tools file_ext
star_genome_index <- function(genome_path = NULL,
                              gff_path = NULL,
                              mappingdir = "results/04_read_mapping",
                              threads = 10) {
    indexdir <- paste0(mapping_dir, "/genomeIndex")
    if(!dir.exists(indexdir)) { dir.create(indexdir, recursive = TRUE) }
    exonparent <- "Parent"
    if(tools::file_ext(gff_path) == "gtf") { exonparent <- "transcript_id" }
    args <- c("--runThreadN", threads, "--runMode genomeGenerate",
              "--genomeDir", indexdir, "--genomeFastaFiles", genome_path,
              "--sjdbGTFfile", gff_path, 
              "--sjdbGTFtagExonParentTranscript", exonparent, 
              "--sjdbOverhang 100")
    system2("STAR", args = args)
    return(NULL)
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
#' @param gff_path Path to the .gff/.gtf file with annotations.
#' @param threads Number of threads for STAR aligner.
#' @importFrom tools file_ext
star_align <- function(sample_info = NULL,
                       filtdir = "results/03_filtered_FASTQ",
                       fastqc_table = NULL,
                       mappingdir = "results/04_read_mapping",
                       gff_path = NULL,
                       threads = 20) {
    exonparent <- "Parent"
    if(tools::file_ext(gff_path) == "gtf") { exonparent <- "transcript_id" }
    fastqc_table$Sample <- gsub("_.*", "", fastqc_table$Sample)
    qc_table <- fastqc_table[!duplicated(fastqc_table$Sample), 
                             c("Sample", "Sequence.length")]
    sample_info2 <- merge(sample_info[!duplicated(sample_info$BioSample),],
                          qc_table, by.x="Run", "Sample")
    indexdir <- paste0(mapping_dir, "/genomeIndex")
    reads <- star_reads(sample_info, filtdir)
    aln <- lapply(seq_len(nrow(sample_info2)), function(x) {
        platform <- sample_info[x, "Instrument"]
        if(grepl("SOLiD|PacBio", platform)) {
            message("Skipping PacBio/SOLiD reads...")
        } else {
            prefix <- paste0(mappingdir, "/", sample_info2[x, "BioSample"])
            seqlen <- sample_info2[x, "Sequence.length"]
            sjdbover <- as.numeric(seqlen) - 1
            args <- c("--runThreadN", threads, "--genomeDir", indexdir, 
                      "--readFilesIn", reads[[x]], "--readFilesCommand zcat",
                      "--outFileNamePrefix", prefix, 
                      "--outSAMtype BAM SortedByCoordinate",
                      "--outSAMstrandField intronMotif", 
                      "--sjdbGTFfile", gff_path, 
                      "--sjdbGTFtagExonParentTranscript", exonparent,
                      "--sjdbOverhang", sjdbover)
            system2("STAR", args = args)
        }
    })
    return(NULL)
}





