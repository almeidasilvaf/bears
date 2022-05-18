
#----Set up conda environemnt management with Herper----------------------------
library(Herper)
my_miniconda <- file.path(tempdir(), "miniconda_bears")

#----Create conda environments for each external dependency
## MultiQC
multiqc_env <- install_CondaTools("multiqc", "multiqc_env", 
                                  pathToMiniConda = my_miniconda)

system2("source", args = "env/miniconda_atlas/etc/profile.d/conda.sh")
system2("conda", args = "init")
system2("source", args = "~./bashrc")

with_CondaEnv("multiqc_env",
              system2("multiqc", args = "-h"),
              pathToMiniConda = my_miniconda)

## FastQC
fastqc_env <- install_CondaTools("fastqc", "fastqc_env", 
                                 pathToMiniConda = my_miniconda)

with_CondaEnv("fastqc_env",
              system2("fastqc", args = "-h"),
              pathToMiniConda = my_miniconda)


## Trimmomatic
trimmomatic_env <- install_CondaTools("trimmomatic", "trimmomatic_env", 
                                      pathToMiniConda = my_miniconda)

## SortMeRNA
sortmerna_env <- install_CondaTools("sortmerna", "sortmerna_env", 
                                    pathToMiniConda = my_miniconda)

## STAR
star_env <- install_CondaTools("star", "star_env", 
                               pathToMiniConda = my_miniconda)

## RSeQC
rseqc_env <- install_CondaTools("rseqc", "rseqc_env", 
                                pathToMiniConda = my_miniconda)

## salmon
salmon_env <- install_CondaTools("salmon", "salmon_env", 
                                 pathToMiniConda = my_miniconda)


## kallisto
kallisto_env <- install_CondaTools("kallisto", "kallisto_env", 
                                   pathToMiniConda = my_miniconda)

## Subread
subread_env <- install_CondaTools("subread", "subread_env", 
                                  pathToMiniConda = my_miniconda)

## StringTie
stringtie_env <- install_CondaTools("stringtie", "stringtie_env", 
                                    pathToMiniConda = my_miniconda)

## TACO
taco_env <- install_CondaTools("taco", "taco_env", 
                               pathToMiniConda = my_miniconda)





# Create conda environment with herper
# envname: "bear_env"
# miniconda_path: my_miniconda


my_miniconda <- file.path("/Users/almeidasilvaf", "miniconda_herper")
Herper::install_CondaTools("multiqc", "bear_env", pathToMiniConda = my_miniconda)
Herper::install_CondaTools("sra-tools==2.10.1", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda) # 2.11.0
Herper::install_CondaTools("fastqc", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("trimmomatic", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("sortmerna", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("star", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("rseqc", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda)
# Herper::install_CondaTools("shrimp", "bear_env", channels = "biobuilds",
#                            updateEnv = TRUE, pathToMiniConda = my_miniconda)
Herper::install_CondaTools("salmon", "bear_env", updateEnv = TRUE,
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("kallisto", "bear_env", updateEnv = TRUE,
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("subread", "bear_env", updateEnv = TRUE,
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("stringtie", "bear_env", updateEnv = TRUE,
                           pathToMiniConda = my_miniconda)

# Create new env for TACO
Herper::install_CondaTools("taco", "taco_env", pathToMiniConda = my_miniconda,
                           updateEnv = TRUE)

