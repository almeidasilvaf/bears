

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
              "--sjdbOverhang 99")
    system2("STAR", args = args)
    return(NULL)
}









