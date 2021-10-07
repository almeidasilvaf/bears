## Code to prepare data in data/
#----sample_info----
term <- "SAMN02422669[BSPL]"
sample_info <- create_sample_info(term)
sample_info$Study_abstract <- textclean::replace_non_ascii(
    sample_info$Study_abstract
)
sample_info$Orientation <- "first"

#----fastqc_table----
data(sample_info)
fq <- system.file("extdata", package="bears")
out <- tempdir()
fqc <- tempdir()
run_fastqc(sample_info, fastqdir = fq, fastqcdir = fqc,
           envname = "bear_env", miniconda_path = my_miniconda)
fastqc_table <- multiqc(fqc, out, envname = "bear_env", 
                        miniconda_path = my_miniconda)

#----mapping_qc----
data(sample_info)
dir <- system.file("extdata", package="bears")
out <- tempdir()
mapping_qc <- multiqc(dir, out, runon="star")


#----tx2gene----
library(GenomicFeatures)
gtf <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf", package="bears")
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

#----Create data sets----
usethis::use_data(sample_info, overwrite = TRUE, compress="xz")
usethis::use_data(fastqc_table, overwrite = TRUE, compress="xz")
usethis::use_data(mapping_qc, overwrite = TRUE, compress="xz")
usethis::use_data(tx2gene, overwrite = TRUE, compress="xz")



###################################################
# Code to prepare data in inst/extdata
#----SRR1039508.fastq and SRR1039508.bam----
library(here)

bampath1 <- system.file("extdata", "SRR1039508_subset.bam", package="airway")
file.copy(from = bampath1, to = here("inst", "extdata"))
# Bash:
# cd inst/extdata
# samtools sort -n SRR1039508.bam -o SRR1039508_sorted.bam
# samtools fastq -@ 8 SRR1039508_sorted.bam \
#   -1 SRR1039508_1.fastq.gz \
#   -2 SRR1039508_2.fastq.gz \
#   -0 /dev/null -s /dev/null -n

#----Homo_sapiens.GRCh37.75_subset.gtf----
gtf_hsapiens <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf",
                            package="airway")
file.copy(from = gtf_hsapiens, to = here("inst", "extdata"), overwrite = TRUE)

#----Homo_sapiens.GRCh37.75_subset.bed----
gtf_hsapiens <- here("inst", "extdata", "Homo_sapiens.GRCh37.75_subset.gtf")
gtf <- rtracklayer::import(gtf_hsapiens)
gtf$score <- as.numeric(rep(0, length(gtf)))
bedfile <- gsub(".gtf", ".bed", gtf_hsapiens)
rtracklayer::export.bed(gtf, bedfile)


#----SRR6967125.fastq.gz----
my_miniconda <- file.path("/Users/almeidasilvaf", "miniconda_herper")
data(sample_info)
fastqdir <- "~/Documents/bear/results/01_FASTQ_files"
sradir <- "~/Documents/bear/results/00_SRA_files"
download_fastq(sample_info[1, ], 
               fastqdir = fastqdir,
               sradir = sradir,
               threads = 2,
               envname = "bear_env", 
               miniconda_path = my_miniconda)
system2("cat", 
        args = c(
            "~/Documents/bear/results/01_FASTQ_files/SRR3725560.fastq.gz",
            " | zcat | head -n 1000",
            "> ~/Documents/bear/results/01_FASTQ_files/SRR3725560.fastq"
            )
        )

#----SRR1039508_1_fastqc.zip----
data(sample_info)
fq <- system.file("extdata", package="bear")
fqc <- tempdir()
run_fastqc(sample_info, fastqdir = fq, fastqcdir = fqc)
file.copy(from = paste0(fqc, "/SRR1039508_1_fastqc.zip"), 
          to = here::here("inst", "extdata"))


#----Hsapiens_GRCh37.75_subset.fa----
gtf_hsapiens <- here("inst", "extdata", "Homo_sapiens.GRCh37.75_subset.gtf")
gtf <- rtracklayer::import(gtf_hsapiens)
start <- min(GenomicRanges::start(gtf)) - 50
end <- max(GenomicRanges::end(gtf)) + 50
chr1 <- Biostrings::readDNAStringSet(
    "~/Downloads/Homo_sapiens.GRCh37.dna_sm.chromosome.1.fa"
    )
chr_subset <- Biostrings::subseq(chr1, start = start, end = end)
Biostrings::writeXStringSet(
    chr_subset, 
    filepath = here::here("inst", "extdata", "Hsapiens_GRCh37.75_subset.fa")
)


#----bac_16s_subset.fa----
bac_16s <- Biostrings::readDNAStringSet("https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-bac-16s-id90.fasta")
bac_16s_subset <- bac_16s[1:50]
Biostrings::writeXStringSet(bac_16s_subset, 
                            filepath = "inst/extdata/bac_16s_subset.fa")


#----SAMN05300278.bam----
# Bash: 
# wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.chromosome.1.fa.gz
# gunzip *.gz
# On server:
data(sample_info)
data(fastqc_table)
genome_path <- list.files(pattern=".fa")
gff_path <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf", 
                        package="bears")
mappingdir <- "results"
indexdir <- "results/index"
filtdir <- system.file("extdata", package="bears")
star_genome_index(genome_path, gff_path, mapping_dir, indexdir, 
                  envname="bear_env", miniconda_path = my_miniconda)
star_align(sample_info, filtdir, fastqc_table, mappingdir, 
           indexdir, gff_path, envname="bear_env", miniconda_path = my_miniconda)


#----Hsapiens_GRCh37.75_subset_transcripts.fa----
library(Biostrings)
library(GenomicFeatures)
library(GenomicRanges)
library(BSgenome)
library(Biostrings)
genome <- Biostrings::readDNAStringSet(
    "~/Downloads/Homo_sapiens.GRCh37.dna_sm.chromosome.1.fa"
)
names(genome) <- 1

# Option 1) All exons of a transcript
transcript_ranges <- gtf[gtf$type %in% c("transcript", "exon")]
transcripts <- split(transcript_ranges, transcript_ranges$transcript_id)

# Option 2) All CDSs of a transcript
#transcript_ranges <- gtf[gtf$type %in% c("transcript", "CDS")]
#transcripts <- split(transcript_ranges, transcript_ranges$transcript_id)

tx_seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = genome,
    transcripts = transcripts
)
writeXStringSet(
    tx_seqs,
    filepath = here::here("inst", "extdata", 
                          "Hsapiens_GRCh37.75_subset_transcripts.fa")
)


#----SAMN02422669----
data(sample_info)
data(fastqc_table)
filtdir <- system.file("extdata", package = "bears")
salmonindex <- tempdir()
salmondir <- tempdir()
transcriptome_path <- system.file(
    "extdata", "Hsapiens_GRCh37.75_subset_transcripts.fa", package="bears"
)
salmon_index(salmonindex, transcriptome_path, 
             envname = "bear_env", miniconda_path = my_miniconda)
salmon_quantify(sample_info, fastqc_table, filtdir, 
                salmonindex, salmondir, envname = "bear_env", 
                miniconda_path = my_miniconda)
fs::dir_copy(file.path(salmondir, "SAMN02422669/"),
             here::here("inst", "extdata"))

filtdir <- system.file("extdata", package = "bears")
kallistoindex <- file.path(tempdir(), "transcripts.idx")
kallistodir <- tempdir()
transcriptome_path <- system.file(
    "extdata", "Hsapiens_GRCh37.75_subset_transcripts.fa", package="bears"
)
kallisto_index(kallistoindex, transcriptome_path, envname = "bear_env",
               miniconda_path = my_miniconda)
kallisto_quantify(sample_info, fastqc_table, filtdir, kallistoindex,
                  kallistodir, envname = "bear_env", 
                  miniconda_path = my_miniconda)
file.copy(file.path(kallistodir, "SAMN02422669", "abundance.tsv"),
          here::here("inst", "extdata", "SAMN02422669"))





