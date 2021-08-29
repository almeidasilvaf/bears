
# Create conda environment with herper
# envname: "bear_env"
# miniconda_path: my_miniconda


my_miniconda <- file.path("/Users/almeidasilvaf", "miniconda_herper")
Herper::install_CondaTools("multiqc", "bear_env", pathToMiniConda = my_miniconda)
Herper::install_CondaTools("sra-tools", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda) # 2.11.0
Herper::install_CondaTools("fastqc", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("trimmomatic", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("sortmerna", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("star", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("star", "bear_env", updateEnv = TRUE, 
                           pathToMiniConda = my_miniconda)
Herper::install_CondaTools("shrimp", "bear_env", channels = "biobuilds",
                           updateEnv = TRUE, pathToMiniConda = my_miniconda)


