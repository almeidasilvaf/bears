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


#----Create data sets----
usethis::use_data(sample_info, overwrite = TRUE, compress="xz")


# Code to prepare data in inst/extdata
#----SRR6967125.fastq.gz----
data(sample_info)
fastqdir <- "~/Documents/bear/results/01_FASTQ_files"
sradir <- "~/Documents/bear/results/00_SRA_files"
download_fastq(sample_info[1, ], 
               fastqdir = fastqdir,
               sradir = sradir,
               threads = 2)
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






