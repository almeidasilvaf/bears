Data acquisition
================

# Data in data/

## sample_info.rda

``` r
term <- "SAMN02422669[BSPL]"
sample_info <- create_sample_info(term)
sample_info$Study_abstract <- textclean::replace_non_ascii(
    sample_info$Study_abstract
)
sample_info$Orientation <- "first"

# Save data
usethis::use_data(sample_info, overwrite = TRUE, compress="xz")
```

## mapping_qc.rda

``` r
data(sample_info)
dir <- system.file("extdata", package="bears")
out <- tempdir()
mapping_qc <- multiqc(dir, out, runon="star")

# Save data
usethis::use_data(mapping_qc, overwrite = TRUE, compress="xz")
```

## tx2gene.rda

``` r
library(GenomicFeatures)
gtf <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf", package="bears")
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Save data
usethis::use_data(tx2gene, overwrite = TRUE, compress="xz")
```

# Data in inst/extdata

## SRR1039508.fastq

``` r
library(here)

bampath1 <- system.file("extdata", "SRR1039508_subset.bam", package="airway")
file.copy(from = bampath1, to = here("inst", "extdata"))
```

``` bash
# Bash:
cd inst/extdata
samtools sort -n SRR1039508.bam -o SRR1039508_sorted.bam
samtools fastq -@ 8 SRR1039508_sorted.bam \
  -1 SRR1039508_1.fastq.gz \
  -2 SRR1039508_2.fastq.gz \
  -0 /dev/null -s /dev/null -n
```

## Homo_sapiens.GRCh37.75_subset.gtf

``` r
gtf_hsapiens <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf",
                            package="airway")

## Include only genes "ENSG00000171819" and "ENSG00000120942" for test reasons
genes <- c("ENSG00000171819", "ENSG00000120942")
gtf <- rtracklayer::import(gtf_hsapiens)
gtf_new <- gtf[gtf$gene_id %in% genes]

## Rescale positions so that gene start is position 1 in the GRanges object
start_pos <- min(start(gtf_new)) - 1
start(gtf_new) <- start(gtf_new) - start_pos
end(gtf_new) <- end(gtf_new) - start_pos

## Export ranges as GTF
rtracklayer::export(
    gtf_new,
    con = here::here("inst", "extdata", "Homo_sapiens.GRCh37.75_subset.gtf"),
    format = "gtf"
)
```

## Homo_sapiens.GRCh37.75_subset.bed

``` r
gtf_hsapiens <- here("inst", "extdata", "Homo_sapiens.GRCh37.75_subset.gtf")
gtf <- rtracklayer::import(gtf_hsapiens)
gtf$score <- as.numeric(rep(0, length(gtf)))
bedfile <- gsub(".gtf", ".bed", gtf_hsapiens)
rtracklayer::export.bed(gtf, bedfile)
```

## Homo_sapiens.GRCh37.75_subset.fa

Here, we will only include the sequences of the genes “ENSG00000171819”
and “ENSG00000120942”, as we did in the example GTF file.

``` r
library(GenomicRanges)
# Load GTF file
gtf_hsapiens <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf",
                            package="airway")
gtf <- rtracklayer::import(gtf_hsapiens) # this contains only chr1
gtf <- gtf[gtf$gene_id %in% c("ENSG00000171819", "ENSG00000120942")]

# Load sequences of chromosome 1
chr1 <- Biostrings::readDNAStringSet(
    "ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.chromosome.1.fa.gz"
)

# Subset regions for the target genes in the GTF file
start_pos <- min(start(gtf))
end_pos <- max(end(gtf)) + 10 # also get +10 nucleotides after gene end

chr_subset <- Biostrings::subseq(chr1, start = start, end = end)
Biostrings::writeXStringSet(
    chr_subset, 
    filepath = here::here("inst", "extdata", 
                          "Homo_sapiens.GRCh37.75_subset.fa")
)
```

## bac_16s_subset.fa

``` r
bac_16s <- Biostrings::readDNAStringSet("https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-bac-16s-id90.fasta")
bac_16s_subset <- bac_16s[1:50]
Biostrings::writeXStringSet(bac_16s_subset, 
                            filepath = "inst/extdata/bac_16s_subset.fa")
```

## Hsapiens_GRCh37.75_subset_transcripts.fa

``` r
library(Biostrings)
library(GenomicFeatures)
library(GenomicRanges)
library(BSgenome)
library(Biostrings)
genome <- Biostrings::readDNAStringSet(
    "ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.chromosome.1.fa.gz"
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
                          "Homo_sapiens.GRCh37.75_subset_transcripts.fa.gz"),
    compress = TRUE
)
```

## SAMN02422669

``` r
data(sample_info)
qc_table <- summary_stats_fastp(system.file("extdata", package = "bears"))
filtdir <- system.file("extdata", package = "bears")
salmonindex <- tempdir()
salmondir <- tempdir()
transcriptome_path <- system.file(
    "extdata", "Hsapiens_GRCh37.75_subset_transcripts.fa", package="bears"
)
salmon_index(salmonindex, transcriptome_path, 
             envname = "bear_env", miniconda_path = my_miniconda)
salmon_quantify(sample_info, qc_table, filtdir, 
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
kallisto_quantify(sample_info, qc_table, filtdir, kallistoindex,
                  kallistodir, envname = "bear_env", 
                  miniconda_path = my_miniconda)
file.copy(file.path(kallistodir, "SAMN02422669", "abundance.tsv"),
          here::here("inst", "extdata", "SAMN02422669"))
```
