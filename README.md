
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bear

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test
coverage](https://codecov.io/gh/almeidasilvaf/bear/branch/main/graph/badge.svg)](https://codecov.io/gh/almeidasilvaf/bear?branch=main)
[![R-CMD-check](https://github.com/almeidasilvaf/bear/workflows/R-CMD-check/badge.svg)](https://github.com/almeidasilvaf/bear/actions)
<!-- badges: end -->

The goal of `bear` is to download RNA-seq data from NCBI SRA, preprocess
them, map to the reference genome, and quantify the expression at the
gene level. The goal of bear is to make RNA-seq data analysis pipelines
reproducible, with a framework built on state-of-the art methods and
softwares.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `bear` from GitHub
with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install("almeidasilvaf/bear")
```

## Code of Conduct

Please note that the `bear` project is released with a [Contributor Code
of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.
