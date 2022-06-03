
library(bears)
library(Herper)
options(timeout = 6000)

#----Setup: create program-specific Conda environment from .yml files-----------
## Install miniconda in a temporary directory
my_miniconda <- file.path(tempdir(), "miniconda")
envs <- list.files(
    system.file("extdata", package = "bears"), pattern = ".yml",
    full.names = TRUE
)

sapply(envs, function(x) {
    import_CondaEnv(x, pathToMiniConda = my_miniconda)
})

## Check if all environments were created
list_CondaEnv(my_miniconda)


#----1) Create directory structure ---------------------------------------------
rootdir <- tempdir() # WARNING: replace with your own directory
ds <- create_dir_structure(rootdir = rootdir)

#----2) Create data frame of sample metadata------------------------------------
term <- "SAMN02422669[BSPL]" # WARNING: replace with your own search term
metadata <- create_sample_info(term)

#----3) Download FASTQ files----------------------------------------------------
download <- download_from_ena(metadata, fastqdir = ds$fastqdir)
download # check status

#----4) QC: Sequence quality and presence of adapters---------------------------
## Run FastQC
with_CondaEnv(
    "fastqc_env",
    fastqc <- run_fastqc(
        metadata, 
        fastqdir = ds$fastqdir, 
        fastqcdir = ds$fastqcdir
    ),
    pathToMiniConda = my_miniconda
)
fastqc # check status

## Run MultiQC
multiqc_out <- file.path(tempdir(), "multiqc")
with_CondaEnv(
    "multiqc_env",
    fastqc_table <- multiqc(ds$fastqcdir, multiqc_out),
    pathToMiniConda = my_miniconda
)
fastqc_table # check summary table

#----5.1) Cleaning 1: Trimming and adapter removal------------------------------
## Run Trimmomatic on reads that failed QC (if any)
with_CondaEnv(
    "trimmomatic_env",
    trim <- trim_reads(metadata, fastqc_table, ds$fastqdir, ds$filtdir),
    pathToMiniConda = my_miniconda
)

#----5.2) Cleaning 2: removal of rRNA-------------------------------------------
## Create rRNA db
rRNA_urls <- c(
    "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/rfam-5.8s-database-id98.fasta",
    "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/rfam-5s-database-id98.fasta",
    "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-arc-16s-id95.fasta",
    "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-arc-23s-id98.fasta",
    "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-bac-16s-id90.fasta",
    "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-bac-23s-id98.fasta",
    "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-euk-18s-id95.fasta",
    "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-euk-28s-id98.fasta"
)
rRNA_db <- file.path(ds$fastqcdir, "rRNA_db")
dir.create(rRNA_db)
lapply(rRNA_urls, function(x) {
    filename <- file.path(rRNA_db, basename(x))
    download.file(x, destfile = filename)
})

## Run SortMeRNA
with_CondaEnv(
    "sortmerna_env",
    sortmerna <- remove_rrna(
        metadata, 
        fastqdir = ds$fastqdir,
        filtdir = ds$filtdir,
        rrna_db_dir = rRNA_db
    ),
    pathToMiniConda = my_miniconda
)

#----6) Read mapping to reference genome----------------------------------------
## Index genome
genome_path <- system.file("extdata", "Hsapiens_GRCh37.75_subset.fa", 
                           package = "bears")
gff_path <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf", 
                        package = "bears")
with_CondaEnv(
    "star_env",
    gen_idx <- star_genome_index(genome_path, gff_path, ds$mapping_dir),
    pathToMiniConda = my_miniconda
)

## Map reads to the genome
with_CondaEnv(
    "star_env",
    mapping <- star_align(
        metadata, filtdir = ds$filtdir, 
        fastqc_table = fastqc_table, 
        mappingdir = mappingdir, 
        gff_path = gff_path
    ),
    pathToMiniconda = my_miniconda
)














