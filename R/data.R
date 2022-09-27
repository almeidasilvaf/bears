
#' Sample metadata
#'
#' This data frame was created with \code{create_sample_info}, and it contains
#' metadata associated with the BioSample SAMN02422669, which is part of the
#' data in the airway package. An additional column named "Orientation" was
#' included with \code{infer_strandedness}, representing read orientation.
#'
#' @name sample_info
#' @format A data frame with sample metadata created 
#' with \code{create_sample_info} and \code{infer_strandedness}.
#' @examples
#' data(sample_info)
#' @usage data(sample_info)
"sample_info"

#' Summary statistics of STAR mapping QC
#' 
#' This data frame was created with \code{multiqc}. The code can be found in 
#' the script datasets.R
#' 
#' @name mapping_qc
#' @format A data frame with read mapping QC summary statistics.
#' @examples 
#' data(mapping_qc)
#' @usage data(mapping_qc)
"mapping_qc"


#' Transcript to gene mapping
#' 
#' The code used to create this data frame can be found in 
#' the script datasets.R
#' 
#' @name tx2gene
#' @format A data frame with transcript IDs in the first column and gene IDs
#' in the second column.
#' @examples 
#' data(tx2gene)
#' @usage data(tx2gene)
"tx2gene"


