
#' Wrapper to get paths to .fastq files depending on the layout
#' 
#' @param fastqdir Path to the directory where .fastq files will be stored.
#' @param run Run accession.
#' @param cmd What command to run. One of "S", "P1", or "P2". S will
#' run the command for single-end reads. P1 will run the command for 
#' the first pair of paired-end reads. P2 will run the command for 
#' the second pair.
#' @return Character with the path to .fastq file.
#' @noRd
get_fastq_paths <- function(fastqdir, run, cmd = "S") {
    
    fastqfile <- NULL
    if(cmd == "S") {
        if(file.exists(paste0(fastqdir, "/", run, ".fastq.gz"))) {
            fastqfile <- paste0(fastqdir, "/", run, ".fastq.gz")
        } else { message("File not found") }
    } else if(cmd == "P1") {
        if(file.exists(paste0(fastqdir, "/", run, "_1.fastq.gz"))) {
            fastqfile <- paste0(fastqdir, "/", run, "_1.fastq.gz")
        } else { message("File not found.") }
    } else {
        if(file.exists(paste0(fastqdir, "/", run, "_2.fastq.gz"))) {
            fastqfile <- paste0(fastqdir, "/", run, "_2.fastq.gz")
        } else { message("File not found.") }
    }
    return(fastqfile)
}

#' Get mapping summary statistics from STAR
#' 
#' @param star_dir Directory where STAR .log files are stored.
#' Default: results/03_read_mapping.
#' 
#' @return A data frame with STAR summary stats for each sample containing the
#' following variables (all numeric, except \strong{Sample}):
#' * Sample
#' * total_reads 
#' * avg_input_read_length
#' * uniquely_mapped
#' * uniquely_mapped_percent
#' * avg_mapped_read_length
#' * num_splices
#' * num_annotated_splices
#' * num_GTAG_splices
#' * num_GCAG_splices
#' * num_ATAC_splices
#' * num_noncanonical_splices
#' * mismatch_rate
#' * deletion_rate
#' * deletion_length
#' * insertion_rate
#' * insertion_length
#' * multimapped
#' * multimapped_percent
#' * multimapped_toomany
#' * multimapped_toomany_percent
#' * unmapped_mismatches_percent
#' * unmapped_tooshort_percent
#' * unmapped_other_percent
#' * unmapped_mismatches
#' * unmapped_tooshort
#' * unmapped_other
#' @export
#' @rdname summary_stats_star
#' @examples
#' star_dir <- system.file("extdata", package = "bears")
#' qc_table <- summary_stats_star(star_dir)
summary_stats_star <- function(star_dir = "results/03_read_mapping") {
    
    star_output <- list.files(
        star_dir, pattern = "Log.final.out$", full.names = TRUE
    )
    
    val <- function(df, p) {
        res <- df$value[grep(p, df$key)]
        if(endsWith(res, "%")) { res <- gsub("%", "", res) }
        return(as.numeric(res))
    }

    parsed <- Reduce(rbind, lapply(seq_along(star_output), function(x) {
        l <- readLines(star_output[x])
        sl <- strsplit(l, " \\|\\\t")
        sl <- split_l[vapply(sl, length, numeric(1)) == 2]
        df <- Reduce(rbind, lapply(sl, function(y) {
            return(data.frame(key = y[1], value = y[2]))
        }))
        
        # Create a data frame of summary stats for sample x
        stats_df <- data.frame(
            Sample = gsub("Log.*", "", basename(star_output[x])),
            total_reads = val(df, "input reads"),
            avg_input_read_length = val(df, "Average input"),
            uniquely_mapped = val(df, "Uniquely mapped reads number"),
            uniquely_mapped_percent = val(df, "Uniquely mapped reads %"),
            avg_mapped_read_length = val(df, "Average mapped length"),
            num_splices = val(df, "Number of splices: Total"),
            num_annotated_splices = val(df, "Number of splices: Annotated"),
            num_GTAG_splices = val(df, "GT/AG"),
            num_GCAG_splices = val(df, "GC/AG"),
            num_ATAC_splices = val(df, "AT/AC"),
            num_noncanonical_splices = val(df, "Non-canonical"),
            mismatch_rate = val(df, "Mismatch rate"),
            deletion_rate = val(df, "Deletion rate"),
            deletion_length = val(df, "Deletion average length"),
            insertion_rate = val(df, "Insertion rate"),
            insertion_length = val(df, "Insertion average length"),
            multimapped = val(df, "Number of reads mapped to multiple"),
            multimapped_percent = val(df, "% of reads mapped to multiple"),
            multimapped_toomany = val(df, "Number of reads mapped to too"),
            multimapped_toomany_percent = val(df, "% of reads mapped to too"),
            unmapped_mismatches_percent = val(df, "% of reads unmapped: too many"),
            unmapped_tooshort_percent = val(df, "% of reads unmapped: too short"),
            unmapped_other_percent = val(df, "% of reads unmapped: other"),
            unmapped_mismatches = val(df, "Number of reads unmapped: too many mis"),
            unmapped_tooshort = val(df, "Number of reads unmapped: too short"),
            unmapped_other = val(df, "Number of reads unmapped: other")
        )
        return(stats_df)
    }))
    return(parsed)
}


#' Get read quality summary statistics from fastp
#'
#' @param fastp_qcdir Character with path to the directory where .json files
#' from fastp are stored. Default: results/QC_dir/fastp_stats.
#'
#' @return A data frame of fastp summary stats for each sample with the 
#' following variables:
#' * Sample
#' * sequencing
#' * before_nreads
#' * before_nbases
#' * before_q20bases
#' * before_q30bases
#' * before_q20rate
#' * before_q30rate
#' * before_GCcontent
#' * before_meanlength
#' * after_nreads
#' * after_nbases
#' * after_q20bases
#' * after_q30bases
#' * after_q20rate
#' * after_q30rate
#' * after_GCcontent
#' * after_meanlength
#' * filter_n_passed
#' * filter_n_lowquality
#' * filter_n_too_many_N
#' * filter_n_tooshort
#' * filter_n_toolong
#' * duplication_rate
#' @importFrom jsonlite fromJSON
#' @export
#' @rdname summary_stats_fastp
#' @examples 
#' fastp_qcdir <- system.file("extdata", package = "bears")
#' fastp_stats <- summary_stats_fastp(fastp_qcdir)
summary_stats_fastp <- function(fastp_qcdir = "results/QC_dir/fastp_stats") {
    
    file <- list.files(fastp_qcdir, pattern = "\\.json$", full.names = TRUE)

    parsed <- Reduce(rbind, lapply(seq_along(file), function(x) {
        j <- jsonlite::fromJSON(file[x])
        js <- as.data.frame(j$summary)
        jf <- as.data.frame(j$filtering_result)
        
        before_length_idx <- intersect(
            grep("mean_length", names(js)), grep("before", names(js))
        )[1]
        after_length_idx <- intersect(
            grep("mean_length", names(js)), grep("after", names(js))
        )[1]
        
        stats_df <- data.frame(
            Sample = gsub("\\.json", "", basename(file[x])),
            sequencing = js$sequencing,
            # Before
            before_nreads = js$before_filtering.total_reads,
            before_nbases = js$before_filtering.total_bases,
            before_q20bases = js$before_filtering.q20_bases,
            before_q30bases = js$before_filtering.q30_bases,
            before_q20rate = js$before_filtering.q20_rate,
            before_q30rate = js$before_filtering.q30_rate,
            before_GCcontent = js$before_filtering.gc_content,
            before_meanlength = js[[before_length_idx]],
            # After
            after_nreads = js$after_filtering.total_reads,
            after_nbases = js$after_filtering.total_bases,
            after_q20bases = js$after_filtering.q20_bases,
            after_q30bases = js$after_filtering.q30_bases,
            after_q20rate = js$after_filtering.q20_rate,
            after_q30rate = js$after_filtering.q30_rate,
            after_GCcontent = js$after_filtering.gc_content, 
            after_meanlength = js[[after_length_idx]],
            # Filtering summary
            filter_n_passed = jf$passed_filter_reads,
            filter_n_lowquality = jf$low_quality_reads,
            filter_n_too_many_N = jf$too_many_N_reads,
            filter_n_tooshort = jf$too_short_reads,
            filter_n_toolong = jf$too_long_reads,
            duplication_rate = j$duplication$rate
        )
        return(stats_df)
    }))
    
    return(parsed)
}




