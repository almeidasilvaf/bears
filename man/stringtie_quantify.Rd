% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05_04_assembly_quantification_stringtie.R
\name{stringtie_quantify}
\alias{stringtie_quantify}
\title{Quantify expression in TPM with StringTie}
\usage{
stringtie_quantify(
  sample_info = NULL,
  qc_table = NULL,
  mappingdir = "results/04_read_mapping",
  gff_path = NULL,
  stringtiedir = "results/05_quantification/stringtie",
  threads = 1
)
}
\arguments{
\item{sample_info}{Data frame of sample metadata created with the
functions \code{create_sample_info} and \code{infer_strandedness}.
The function \code{infer_strandedness} adds a column named "Orientation"
with library strandedness information, which is mandatory.}

\item{qc_table}{Data frame of fastp summary statistics as returned
by \code{summary_stats_fastp()}.}

\item{mappingdir}{Directory where .bam files are stored.}

\item{gff_path}{Path to GFF/GTF file with annotations.}

\item{stringtiedir}{Directory where StringTie output files will be stored.}

\item{threads}{Number of threads to use. Default: 1.}
}
\value{
A 2-column data frame with samples in the first column and status
in the second column, with "OK" if expression in TPM was obtained, and
NA otherwise.
}
\description{
Quantify expression in TPM with StringTie
}
\examples{
data(sample_info)
qc_table <- summary_stats_fastp(system.file("extdata", package = "bears"))
mappingdir <- system.file("extdata", package="bears")
gff_path <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf",
                        package="bears")
stringtiedir <- tempdir()
if(stringtie_is_installed()) {
    a <- stringtie_quantify(sample_info, qc_table, mappingdir, 
                            gff_path, stringtiedir)
}
}
