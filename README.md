
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GEAtlas

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check-bioc](https://github.com/almeidasilvaf/GEAtlas/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/almeidasilvaf/GEAtlas/actions)
[![Codecov test
coverage](https://codecov.io/gh/almeidasilvaf/GEAtlas/branch/main/graph/badge.svg)](https://codecov.io/gh/almeidasilvaf/GEAtlas?branch=main)
<!-- badges: end -->

The goal of `GEAtlas` is to download RNA-seq data from NCBI SRA,
preprocess them, map to the reference genome, and quantify the
expression at the gene level. The goal of GEAtlas is to make RNA-seq
data analysis pipelines reproducible, with a framework built on
state-of-the art methods and softwares.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `GEAtlas` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("GEAtlas")
```

And the development version from
[GitHub](https://github.com/almeidasilvaf/GEAtlas) with:

``` r
remotes::install_github("almeidasilvaf/GEAtlas")
```

## Code of Conduct

Please note that the `GEAtlas` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

-   Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.13/BiocCheck)*.
-   Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation website](http://almeidasilvaf.github.io/GEAtlas)
    is automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
-   The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.13/biocthis)*.
