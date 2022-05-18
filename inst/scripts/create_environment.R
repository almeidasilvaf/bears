
#----Set up conda environemnt management with Herper----------------------------
library(Herper)
library(here)
my_miniconda <- file.path(tempdir(), "miniconda_bears")

#----Create conda environments for each external dependency
## MultiQC
multiqc_env <- install_CondaTools("multiqc", "multiqc_env", 
                                  pathToMiniConda = my_miniconda)

export_CondaEnv("multiqc_env",
                yml_export = here("inst", "extdata", "multiqc_env.yml"),
                my_miniconda)

# Test that it works
# with_CondaEnv("multiqc_env",
#               system2("multiqc", args = "-h"),
#               pathToMiniConda = my_miniconda)

## FastQC
fastqc_env <- install_CondaTools("fastqc", "fastqc_env", 
                                 pathToMiniConda = my_miniconda)

export_CondaEnv("fastqc_env",
                yml_export = here("inst", "extdata", "fastqc_env.yml"),
                my_miniconda)

## Trimmomatic
trimmomatic_env <- install_CondaTools("trimmomatic", "trimmomatic_env", 
                                      pathToMiniConda = my_miniconda)

export_CondaEnv("trimmomatic_env",
                yml_export = here("inst", "extdata", "trimmomatic_env.yml"),
                my_miniconda)

## SortMeRNA
sortmerna_env <- install_CondaTools("sortmerna", "sortmerna_env", 
                                    pathToMiniConda = my_miniconda)

export_CondaEnv("sortmerna_env",
                yml_export = here("inst", "extdata", "sortmerna_env.yml"),
                my_miniconda)

## STAR
star_env <- install_CondaTools("star", "star_env", 
                               pathToMiniConda = my_miniconda)

export_CondaEnv("star_env",
                yml_export = here("inst", "extdata", "star_env.yml"),
                my_miniconda)

## RSeQC
rseqc_env <- install_CondaTools("rseqc", "rseqc_env", 
                                pathToMiniConda = my_miniconda)

export_CondaEnv("rseqc_env",
                yml_export = here("inst", "extdata", "rseqc_env.yml"),
                my_miniconda)

## salmon
salmon_env <- install_CondaTools("salmon", "salmon_env", 
                                 pathToMiniConda = my_miniconda)

export_CondaEnv("salmon_env",
                yml_export = here("inst", "extdata", "salmon_env.yml"),
                my_miniconda)

## kallisto
kallisto_env <- install_CondaTools("kallisto", "kallisto_env", 
                                   pathToMiniConda = my_miniconda)

export_CondaEnv("kallisto_env",
                yml_export = here("inst", "extdata", "kallisto_env.yml"),
                my_miniconda)

## Subread
subread_env <- install_CondaTools("subread", "subread_env", 
                                  pathToMiniConda = my_miniconda)

export_CondaEnv("subread_env",
                yml_export = here("inst", "extdata", "subread_env.yml"),
                my_miniconda)


## StringTie
stringtie_env <- install_CondaTools("stringtie", "stringtie_env", 
                                    pathToMiniConda = my_miniconda)

export_CondaEnv("stringtie_env",
                yml_export = here("inst", "extdata", "stringtie_env.yml"),
                my_miniconda)

## TACO
taco_env <- install_CondaTools("taco", "taco_env", 
                               pathToMiniConda = my_miniconda)

export_CondaEnv("taco_env",
                yml_export = here("inst", "extdata", "taco_env.yml"),
                my_miniconda)
