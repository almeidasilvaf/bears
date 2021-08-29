## Code to prepare data in data/
#----sample_info----
term_se <- "PRJNA324514[BPRJ]"
df_se <- create_sample_info(term_se)
df_se <- df_se[2, ]

# Paired-end example data
term_pe <- "SAMD00117551[BSPL]"
df_pe <- create_sample_info(term_pe)

# SOLiD example data
term_solid <- "SAMN02688407[BSPL]"
df_solid <- create_sample_info(term_solid)

# Combine data frames
sample_info <- rbind(df_se, df_pe, df_solid)
sample_info$Tissue <- c("leaf", "leaf", "root")

#----fastqc_table----
data(sample_info)
fq <- system.file("extdata", package="bear")
out <- tempdir()
fastqc_table <- multiqc(fq, out, envname = "bear_env", 
                        miniconda_path = my_miniconda)

#----Create data sets----
usethis::use_data(sample_info, overwrite = TRUE, compress="xz")
usethis::use_data(fastqc_table, overwrite = TRUE, compress="xz")


###################################################
# Code to prepare data in inst/extdata
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

#----SRR3725560_fastqc.zip----
data(sample_info)
fq <- system.file("extdata", package="bear")
fqc <- tempdir()
run_fastqc(sample_info[1, ], fastqdir = fq, fastqcdir = fqc)


#----bac_16s_subset.fa----
bac_16s <- Biostrings::readDNAStringSet("https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-bac-16s-id90.fasta")
bac_16s_subset <- bac_16s[1:50]
Biostrings::writeXStringSet(bac_16s_subset, 
                            filepath = "inst/extdata/bac_16s_subset.fa")

#----Gmax_chr15_subset.fa----
Gmax_chr15 <- Biostrings::readDNAStringSet("ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/Genomes/gma.con.gz")
Gmax_chr15_subset <- Gmax_chr15[names(Gmax_chr15) == "Chr15"]
Gmax_chr15_subset2 <- Biostrings::subseq(Gmax_chr15_subset, start=1, end = 100000)
Biostrings::writeXStringSet(Gmax_chr15_subset2,
                            filepath = "inst/extdata/Gmax_chr15_subset.fa")


#----Gmax_chr15_subset.gff3----
Gmax_chr15_ranges <- rtracklayer::import("ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/GFF/gma/annotation.all_transcripts.all_features.gma.gff3.gz")
Gmax_chr15_ranges_subset <- Gmax_chr15_ranges[seqnames(Gmax_chr15_ranges) == "Chr15"]
Gmax_chr15_ranges_subset <- Gmax_chr15_ranges_subset[1:268]
rtracklayer::export.gff3(Gmax_chr15_ranges_subset,
                         con = "inst/extdata/Gmax_chr15_subset.gff3")


#----SAMN05300278.bam----
data(sample_info)
data(fastqc_table)
genome_path <- system.file("extdata", "Gmax_chr15_subset.fa", package="bears")
gff_path <- system.file("extdata", "Gmax_chr15_subset.gff3", package="bears")
mappingdir <- "results/"
indexdir <- "results/index"
filtdir <- system.file("extdata", package="bears")
star_genome_index(genome_path, gff_path, mapping_dir, indexdir, 
                  envname="bear_env", miniconda_path = my_miniconda)
star_align(sample_info[1, ], filtdir, fastqc_table, mappingdir, 
           indexdir, gff_path, envname="bear_env", miniconda_path = my_miniconda)


