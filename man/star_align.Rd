% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_read_mapping.R
\name{star_align}
\alias{star_align}
\title{Align reads to a reference genome using STAR}
\usage{
star_align(
  sample_info = NULL,
  filtdir = "results/03_filtered_FASTQ",
  fastqc_table = NULL,
  mappingdir = "results/04_read_mapping",
  indexdir = "results/04_read_mapping/genomeIndex",
  gff_path = NULL,
  threads = 20
)
}
\arguments{
\item{sample_info}{Data frame of sample metadata created with the
function \code{create_sample_info}.}

\item{filtdir}{Path to the directory where filtered reads will be stored.
Default: results/03_filtered_FASTQ.}

\item{fastqc_table}{Data frame of summary statistics for FastQC as returned
by \code{multiqc()}.}

\item{mappingdir}{Path to the directory where read mapping files (.bam) will
be stored.}

\item{indexdir}{Directory where the STAR genome index files will be stored.
Default: results/04_read_mapping/genomeIndex.}

\item{gff_path}{Path to the .gff/.gtf file with annotations.}

\item{threads}{Number of threads for STAR aligner.}
}
\description{
Align reads to a reference genome using STAR
}