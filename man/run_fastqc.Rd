% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_qc.R
\name{run_fastqc}
\alias{run_fastqc}
\title{Run FastQC on .fastq files}
\usage{
run_fastqc(
  sample_info,
  fastqdir = "results/01_FASTQ_files",
  fastqcdir = "results/02_FastQC_dir"
)
}
\arguments{
\item{sample_info}{Data frame of sample metadata created with the
function \code{create_sample_info}.}

\item{fastqdir}{Path to the directory where .fastq files will be stored.
Default: results/01_FASTQ_files.}

\item{fastqcdir}{Path to the directory where FastQC output will be stored.
Default: results/02_FastQC_dir.}
}
\value{
A 2-column data frame with run accession in the first column and
FastQC run status. If FastQC ran successfully, the status "OK" is displayed.
If FastQC failed to process a particular file, NA will be displayed on its
corresponding name.
}
\description{
Run FastQC on .fastq files
}
\examples{
data(sample_info)
fq <- system.file("extdata", package="bears")
fqc <- tempdir()
if(fastqc_is_installed()) {
    run_fastqc(sample_info, fastqdir = fq, fastqcdir = fqc)
}
}
