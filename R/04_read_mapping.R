
#' Index genome for STAR alignment
#'
#' @param genome_path Path to the genome .fasta file.
#' @param gff_path Path to the .gff/.gtf file with annotations.
#' @param mappingdir Path to the directory where read mapping files (.bam) will
#' be stored.
#' @param indexdir Directory where the STAR genome index files will be stored.
#' Default: results/04_read_mapping/genomeIndex.
#' @param threads Number of threads for STAR aligner.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#' 
#' @return NULL
#' @export
#' @rdname star_genome_index
#' @importFrom tools file_ext
#' @examples 
#' \donttest{
#' genome_path <- system.file("extdata", "Gmax_chr15_subset.fa", package="bear")
#' gff_path <- system.file("extdata", "Gmax_chr15_subset.gff3", package="bear")
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
                              threads = 10,
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
#' @param indexdir Directory where the STAR genome index files will be stored.
#' Default: results/04_read_mapping/genomeIndex.
#' @param gff_path Path to the .gff/.gtf file with annotations.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#'
#' @return NULL
#' @param threads Number of threads for STAR aligner.
#' @importFrom tools file_ext
#' @export
#' @rdname star_align
#' @examples
#' \donttest{
#' data(sample_info)
#' data(fastqc_table)
#' genome_path <- system.file("extdata", "Gmax_chr15_subset.fa", package="bear")
#' gff_path <- system.file("extdata", "Gmax_chr15_subset.gff3", package="bear")
#' mappingdir <- tempdir()
#' indexdir <- tempdir()
#' threads <- 2
#' filtdir <- system.file("extdata", package="bear")
#' if(star_is_installed()) {
#'     star_genome_index(genome_path, gff_path, mapping_dir, indexdir)
#'     star_align(sample_info[1, ], filtdir, fastqc_table, mappingdir, 
#'                indexdir, gff_path, threads)
#' }
#' }
star_align <- function(sample_info = NULL,
                       filtdir = "results/03_filtered_FASTQ",
                       fastqc_table = NULL,
                       mappingdir = "results/04_read_mapping",
                       indexdir = "results/04_read_mapping/genomeIndex",
                       gff_path = NULL,
                       threads = 20,
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
#' @noRd
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



#' Align SOLiD reads with gmapper-cs
#'
#' @param sample_info Data frame of sample metadata.
#' @param genome_path Path to the FASTA file representing the genome.
#' @param soliddir Path to the directory where SOLiD reads are stored.
#' @param mappingdir Path to the directory where .bam files will be stored.
#' @param indexdir Path to the directory where genome index will be stored.
#' @param threads Numeric representing the number of threads to use.
#' Default = 2.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#'
#' @return NULL
#' @importFrom fs file_delete
#' @export
#' @rdname solid_align
#' @examples
#' \donttest{
#' data(sample_info)
#' genome_path <- system.file("extdata", "Gmax_chr15_subset.fa", package="bear")
#' soliddir <- system.file("extdata", package="bear")
#' mappingdir <- tempdir()
#' indexdir <- tempdir()
#' if(shrimp_is_installed()) {
#' solid_align(sample_info[1, ], genome_path, soliddir, mappingdir,
#'             indexdir)
#' }
#' }
solid_align <- function(sample_info = NULL,
                        genome_path = NULL,
                        soliddir = "results/01_SOLiD_dir",
                        mappingdir = "results/04_read_mapping",
                        indexdir = "results/04_read_mapping/shrimpIndex",
                        threads = 2,
                        envname = NULL,
                        miniconda_path = NULL) {
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!shrimp_is_installed()) { stop("Unable to find SHRiMP in PATH.") }
    if(!dir.exists(indexdir)) { dir.create(indexdir, recursive = TRUE) }
    # Create genome index
    args1 <- c("--threads", threads, "--save", indexdir, genome_path)
    system2("gmapper-cs", args = args1)
    # Map reads
    sample_info <- sample_info[grep("SOLiD", sample_info$Instrument), ]
    m <- lapply(seq_len(nrow(sample_info)), function(x) {
        var <- var2list(sample_info, index = x)
        read <- paste0(soliddir, "/", var$run, "*.csfasta")
        outsam <- paste0(mappingdir, "/", var$biosample, "_shrimp.sam")
        outbam <- paste0(mappingdir, "/", var$biosample, "_shrimp.bam")
        outbamsort <- paste0(mappingdir, "/", var$biosample, 
                             "Aligned.sortedByCoord.out.bam")
        margs <- c("--threads", threads, "--local --sam --all-contigs",
                   "--strata --load", indexdir, read, ">", outsam)
        system2("gmapper-cs", margs)
        # SAM to BAM conversion
        sam2bamargs <- c("view -bh -@ 10 -o", outbam, outsam)
        system2("samtools", args = sam2bamargs)
        # BAM sorting
        argssorting <- c("sort -@ 10 -o", outbamsort, outbam)
        system2("samtools", args = argssorting)
        del <- fs::file_delete(c(outbam, outsam))
    })
    return(NULL)
}







