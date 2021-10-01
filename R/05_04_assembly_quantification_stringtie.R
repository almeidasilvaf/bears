

#' Assemble transcripts with StringTie
#' 
#' @param sample_info Data frame of sample metadata created with the
#' functions \code{create_sample_info} and \code{infer_strandedness}.
#' The function \code{infer_strandedness} adds a column named "Orientation" 
#' with library strandedness information, which is mandatory.
#' @param fastqc_table Data frame of summary statistics for FastQC as returned
#' by \code{multiqc()}.
#' @param mappingdir Directory where .bam files are stored.
#' @param gff_path Path to GFF/GTF file with annotations.
#' @param stringtiedir Directory where StringTie output files will be stored.
#' @param threads Number of threads to use. Default: 2.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#' 
#' @return A 2-column data frame with samples in the first column and status
#' in the second column, with "OK" if transcripts were assembled, and 
#' NA otherwise.
#' @importFrom Herper local_CondaEnv
#' @export
#' @rdname stringtie_assemble
#' @examples
#' data(sample_info)
#' data(fastqc_table)
#' mappingdir <- system.file("extdata", package="bears")
#' gff_path <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf",
#'                         package="bears")
#' stringtiedir <- tempdir()
#' if(stringtie_is_installed()) {
#'     a <- stringtie_assemble(sample_info, fastqc_table, mappingdir, 
#'                             gff_path, stringtiedir)
#' }
stringtie_assemble <- function(
    sample_info = NULL, fastqc_table = NULL,
    mappingdir = "results/04_read_mapping",
    gff_path = NULL,
    stringtiedir = "results/05_quantification/stringtie",
    threads = 2,
    envname = NULL, miniconda_path = NULL
    ) {
    
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!stringtie_is_installed()) { stop("Unable to find StringTie in PATH.") }
    if(!dir.exists(stringtiedir)) { dir.create(stringtiedir, recursive = TRUE) }
    
    bamfiles <- list.files(mappingdir, pattern = ".bam", full.names = FALSE)
    bamfiles <- gsub("Aligned.*", "", bamfiles)
    sample_meta <- sample_info[!duplicated(sample_info$BioSample), ]
    sample_meta <- sample_meta[sample_meta$BioSample %in% bamfiles, ]
    fastqc_table$Sample <- gsub("_1", "", fastqc_table$Sample)
    fastqc_table <- fastqc_table[!grepl("_2", fastqc_table$Sample), ]
    
    t <- lapply(seq_len(nrow(sample_meta)), function(x) {
        var <- var2list(sample_meta, index = x)
        infile <- paste0(mappingdir, "/", 
                         var$biosample, "Aligned.sortedByCoord.out.bam")
        outfile <- paste0(stringtiedir, "/assembly/", var$biosample, ".gff")
        orientation <- sample_meta[x, "Orientation"]
        lib <- translate_strandedness(orientation, var$layout)$stringtie
        read_len <- fastqc_table[fastqc_table$Sample == var$biosample, 
                                 "Sequence.length"]
        anchor_len <- read_len / 4
        args <- c(infile, "-G", gff_path, "-o", outfile, 
                  "-p", threads, "-j 5 -c 10 -a", anchor_len, lib)
        system2("stringtie", args = args)
    })
    gff_files <- list.files(paste0(stringtiedir, "/assembly"), full.names=FALSE)
    gff_files <- gsub(".gff", "", gff_files)
    df_status <- data.frame(sample = sample_meta$BioSample)
    df_status$status <- ifelse(df_status$sample %in% gff_files, "OK", NA)
    return(df_status)
}


#' Clean directories after running TACO
#' 
#' @param out2dir Directory where final assembly.gtf file is.
#' @param final_dir Directory where final assembly.gtf file will be stored.
#' 
#' @return A NULL object.
#' @importFrom fs file_move dir_delete file_delete
#' @noRd
taco_clean <- function(out2dir = NULL, final_dir = NULL) {
    
    move <- fs::file_move(paste0(out2dir, "/assembly.gtf"),
                          paste0(final_dir, "/final_assembly.gtf"))
    dir_list <- list.dirs(final_dir, recursive = FALSE, full.names=TRUE)
    dir_list <- dir_list[grepl("biosample_", dir_list)]
    del <- fs::dir_delete(dir_list)
    del <- fs::dir_delete(out2dir)
    txtfiles <- list.files(final_dir, pattern = ".txt", full.names = TRUE)
    txtdel <- fs::file_delete(txtfiles)
    return(NULL)
}

#' Merge assembled transcripts with TACO
#'
#' This function merges GFF files from StringTie assembly into a single GFF.
#' 
#' @param sample_info Data frame of sample metadata created with the
#' function \code{create_sample_info}.
#' @param stringtiedir Directory where StringTie output files will be stored.
#' @param threads Number of threads to use. Default: 2.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#'
#' @details The GFF files for each BioSample are first merged into 
#' BioProject-level GFF files. Then, these Bioproject-level GFF files are
#' merged into a single file that represents the whole set.
#'
#' @return A 2-column data frame with the TACO run status, with "OK" if the
#' program merged transcripts successfully, and NA otherwise.
#' @importFrom Herper local_CondaEnv
#' @export
#' @rdname taco_merge
#' @examples
#' data(sample_info)
#' data(fastqc_table)
#' mappingdir <- system.file("extdata", package="bears")
#' gff_path <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf",
#'                         package="bears")
#' stringtiedir <- tempdir()
#' if(stringtie_is_installed()) {
#'     a <- stringtie_assemble(sample_info, fastqc_table, mappingdir, 
#'                             gff_path, stringtiedir)
#' }
#' if(taco_is_installed()) {
#'     taco_merge(sample_info, stringtiedir)
#' }
taco_merge <- function(
    sample_info = NULL, 
    stringtiedir = "results/05_quantification/stringtie",
    threads = 2,
    envname = NULL, miniconda_path = NULL
    ) {
    
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!taco_is_installed()) { stop("Unable to find TACO in PATH.") }
    assemblydir <- paste0(stringtiedir, "/assembly")
    mergeddir <- paste0(assemblydir, "/merged_assembly")
    if(!dir.exists(mergeddir)) { dir.create(mergeddir, recursive = TRUE ) }
    
    gffs <- list.files(assemblydir, pattern = ".gff", full.names = FALSE)
    gffs <- gsub(".gff", "", gffs)
    sample_meta <- sample_info[!duplicated(sample_info$BioSample), ]
    sample_meta <- sample_meta[sample_meta$BioSample %in% gffs, ]
    
    # Create a .gtf file for each bioproject
    bioproject_list <- split(sample_meta$BioSample, sample_meta$BioProject)
    file_list <- lapply(seq_along(bioproject_list), function(x) {
        files <- paste0(assemblydir, "/", bioproject_list[[x]], ".gff")
        txt_file1 <- paste0(mergeddir, "/", names(bioproject_list)[x], ".txt")
        writeLines(files, txt_file1)
        outdir1 <- paste0(mergeddir, "/biosample_", names(bioproject_list)[x])
        args1 <- c("--num-processes", threads, "-o", outdir1, 
                   "--gtf-expr-attr TPM", txt_file1)
        system2("taco_run", args = args1)
    })
    # Merge everything into a single .gtf file
    bp_gffs <- paste0(mergeddir, "/biosample_", names(bioproject_list),
                      "/assembly.gtf")
    txt_file2 <- paste0(mergeddir, "/all_bioprojects.txt")
    writeLines(bp_gffs, txt_file2)
    outdir2 <- paste0(mergeddir, "/final_assembly")
    args2 <- c("--num-processes", threads, "-o", outdir2, 
               "--gtf-expr-attr expr", txt_file2)
    system2("taco_run", args = args2)
    
    # Cleaning up
    c <- taco_clean(outdir2, mergeddir)
    df_status <- data.frame(file = "final_assembly.gtf",
                            status = NA)
    if(file.exists(paste0(mergeddir, "/final_assembly.gtf"))) {
        df_status$status <- "OK"
    }
    return(df_status)
}


#' Quantify expression in TPM with StringTie
#' 
#' @param sample_info Data frame of sample metadata created with the
#' functions \code{create_sample_info} and \code{infer_strandedness}.
#' The function \code{infer_strandedness} adds a column named "Orientation" 
#' with library strandedness information, which is mandatory.
#' @param fastqc_table Data frame of summary statistics for FastQC as returned
#' by \code{multiqc()}.
#' @param mappingdir Directory where .bam files are stored.
#' @param gff_path Path to GFF/GTF file with annotations.
#' @param stringtiedir Directory where StringTie output files will be stored.
#' @param threads Number of threads to use. Default: 2.
#' @param envname Name of the Conda environment with external dependencies 
#' to be included in the temporary R environment.
#' @param miniconda_path Path to miniconda. Only valid if envname is specified.
#' 
#' @return A 2-column data frame with samples in the first column and status
#' in the second column, with "OK" if expression in TPM was obtained, and 
#' NA otherwise.
#' @importFrom Herper local_CondaEnv
#' @export
#' @rdname stringtie_quantify
#' @examples
#' data(sample_info)
#' data(fastqc_table)
#' mappingdir <- system.file("extdata", package="bears")
#' gff_path <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf",
#'                         package="bears")
#' stringtiedir <- tempdir()
#' if(stringtie_is_installed()) {
#'     a <- stringtie_quantify(sample_info, fastqc_table, mappingdir, 
#'                             gff_path, stringtiedir)
#' }
stringtie_quantify <- function(
    sample_info = NULL, fastqc_table = NULL,
    mappingdir = "results/04_read_mapping",
    gff_path = NULL,
    stringtiedir = "results/05_quantification/stringtie",
    threads = 2,
    envname = NULL, miniconda_path = NULL
) {
    
    if(load_env(envname, miniconda_path)) {
        Herper::local_CondaEnv(envname, pathToMiniConda = miniconda_path)
    }
    if(!stringtie_is_installed()) { stop("Unable to find StringTie in PATH.") }
    if(!dir.exists(stringtiedir)) { dir.create(stringtiedir, recursive = TRUE) }
    quantdir <- paste0(stringtiedir, "/quant/")
    if(!dir.exists(quantdir)) { dir.create(quantdir, recursive = TRUE) }
    
    bamfiles <- list.files(mappingdir, pattern = ".bam", full.names = FALSE)
    bamfiles <- gsub("Aligned.*", "", bamfiles)
    sample_meta <- sample_info[!duplicated(sample_info$BioSample), ]
    sample_meta <- sample_meta[sample_meta$BioSample %in% bamfiles, ]
    fastqc_table$Sample <- gsub("_1", "", fastqc_table$Sample)
    fastqc_table <- fastqc_table[!grepl("_2", fastqc_table$Sample), ]
    
    t <- lapply(seq_len(nrow(sample_meta)), function(x) {
        var <- var2list(sample_meta, index = x)
        infile <- paste0(mappingdir, "/", 
                         var$biosample, "Aligned.sortedByCoord.out.bam")
        outgtf <- paste0(quantdir, var$biosample, ".gtf")
        outbg <- paste0(quantdir, var$biosample, ".bw")
        outabund <- paste0(quantdir, var$biosample, ".gene_abund.out")
        outcov <- paste0(quantdir, var$biosample, ".cov")
        orientation <- sample_meta[x, "Orientation"]
        lib <- translate_strandedness(orientation, var$layout)$stringtie
        read_len <- fastqc_table[fastqc_table$Sample == var$biosample, 
                                 "Sequence.length"]
        anchor_len <- read_len / 4
        args <- c(infile, "-G", gff_path, "-o", outgtf, "-p", threads, 
                  "-e -b", outbg, "-A", outabund, "-C", outcov, 
                  "-j 5 -c 10 -a", anchor_len, lib)
        system2("stringtie", args = args)
    })
    ab_files <- list.files(paste0(stringtiedir, "/quant"), pattern=".out",
                            full.names=FALSE)
    ab_files <- gsub(".gene_abund.out", "", ab_files)
    df_status <- data.frame(sample = sample_meta$BioSample)
    df_status$status <- ifelse(df_status$sample %in% ab_files, "OK", NA)
    return(df_status)
}


#' Create a SummarizedExperiment object from featureCounts output
#' 
#' @param sample_info Data frame of sample metadata created with the
#' functions \code{create_sample_info}
#' @param stringtiedir Directory where StringTie output files are stored.
#' @param level Character indicating to which level expression must be 
#' quantified in the SE object. One of "gene" (default), "transcript", 
#' or "both". For "both", the SE object will have two assays named "transcript"
#' and "gene".
#' @param tx2gene Data frame of correspondence between genes and transcripts, 
#' with gene IDs in the first column and transcript IDs in the second column.
#' Only required if level = 'gene' or 'both'. 
#' 
#' @return A SummarizedExperiment object with gene/transcript expression
#' levels and sample metadata.
#' @importFrom tximport tximport summarizeToGene
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @rdname stringtie2se
#' @examples 
#' data(sample_info)
#' data(fastqc_table)
#' data(tx2gene)
#' mappingdir <- system.file("extdata", package="bears")
#' gff_path <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf",
#'                         package="bears")
#' stringtiedir <- tempdir()
#' if(stringtie_is_installed()) {
#'     a <- stringtie_quantify(sample_info, fastqc_table, mappingdir, 
#'                             gff_path, stringtiedir)
#'     se_gene <- stringtie2se(sample_info, stringtiedir, tx2gene = tx2gene)
#' }
stringtie2se <- function(sample_info = NULL,
                         stringtiedir = "results/05_quantification/stringtie",
                         level = "gene",
                         tx2gene = NULL) {
    quantdir <- paste0(stringtiedir, "/quant")
    dirs <- list.dirs(quantdir, full.names = TRUE, recursive=FALSE)
    dirs <- dirs[grepl("SAMN", dirs)]
    samples <- gsub(paste0(quantdir, "/"), "", dirs)
    samples <- gsub("\\.bw", "", samples)
    files <- paste0(dirs, "/t_data.ctab")
    names(files) <- samples
    sample_meta <- sample_info[!duplicated(sample_info$BioSample), ]
    sample_meta <- sample_meta[sample_meta$BioSample %in% samples, ]
    coldata <- data.frame(row.names = sample_meta$BioSample)
    coldata <- cbind(coldata, sample_meta[, !names(sample_meta) == "BioSample"])
    
    if(level == "gene") {
        exp <- tximport::tximport(files, type = "stringtie", tx2gene = tx2gene)
        final <- SummarizedExperiment::SummarizedExperiment(
            assays = list(gene_TPM = exp$abundance),
            colData = coldata
        )
    } else if(level == "transcript") {
        exp <- tximport::tximport(files, type = "stringtie", txOut = TRUE)
        final <- SummarizedExperiment::SummarizedExperiment(
            assays = list(tx_TPM = exp$abundance),
            colData = coldata
        )
    } else if(level == "both") {
        exp_tx <- tximport::tximport(files, type = "stringtie", txOut = TRUE)
        exp_gene <- tximport::summarizeToGene(exp_tx, tx2gene)
        se_gene <- SummarizedExperiment::SummarizedExperiment(
            assays = list(gene_TPM = exp_gene$abundance),
            colData = coldata
        )
        se_tx <- SummarizedExperiment::SummarizedExperiment(
            assays = list(tx_TPM = exp_tx$abundance),
            colData = coldata
        )
        final <- list(gene = se_gene, transcript = se_tx)
    } else {
        stop("Invalid parameter for the 'level' argument.")
    }
    return(final)
}

