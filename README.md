
<!-- README.md is generated from README.Rmd. Please edit that file -->
SeuratData v0.2.1
=================

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/SeuratData)](https://CRAN.R-project.org/package=SeuratData) [![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://github.com/satijalab/seurat-data) <!-- badges: end -->

SeuratData is a mechanism for distributing datasets in the form of [Seurat](https://satijalab.org/seurat) objects using R's internal package and data management systems. It represents an easy way for users to get access to datasets that are used in the Seurat vignettes.

Installation
------------

SeuratData is not currently available on CRAN. You can install it from GitHub with:

``` r
if (!requireNamespace("remotes", quietly = FALSE)) {
  install.packages("remotes")
}
remotes::install_github("satijalab/seurat-disk")
```

Getting Started
---------------

When loading SeuratData, a list of all available datasets will be displayed (this is similar to other metapackages like [tidyverse](https://cran.r-project.org/package=tidyverse) along with the version of [Seurat](https://satijalab.org/seurat/) used to create each dataset. This message can be suppressed with [`suppressPackageStartupMessages`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/message.html)

``` r
library(SeuratData)
#> Attaching Seurat
#> ── Installed datasets ─────────────────────────────────────────────────────────────── SeuratData v0.2.1 ──
#> x cbmc     3.0.0                                     x stxBrain 0.1.1
#> x pbmc3k   3.0.0
#> ─────────────────────────────────────────────────── Key ──────────────────────────────────────────────────
#> ✓ Dataset built with compatible version of Seurat
#> > Dataset built with a newer version of Seurat than installed
#> x Unknown version of Seurat dataset was built with
```

To see a manifest of all available datasets, use `AvailableData`; this manifest will update as new datasets are uploaded to our data repository.

``` r
AvailableData()
#>          Version SeuratBuild DefaultDS
#> cbmc       3.0.0        <NA>      cbmc
#> hcabm40k   3.0.0        <NA>      raw*
#> ifnb       3.0.0        <NA>      raw*
#> panc8      3.0.2        <NA>      raw*
#> pbmc3k     3.0.0        <NA>    pbmc3k
#> pbmcsca    3.0.0        <NA>      raw*
#> stxBrain   0.1.1        <NA>      <NA>
#>                                               OtherDS DiskDS
#> cbmc                                             <NA>   <NA>
#> hcabm40k                                          NA*    NA*
#> ifnb                                              NA*    NA*
#> panc8                                             NA*    NA*
#> pbmc3k                                   pbmc3k.final   <NA>
#> pbmcsca                                           NA*    NA*
#> stxBrain anterior1, anterior2, posterior1, posterior2   <NA>
#>                                                                 Summary species
#> cbmc                       scRNAseq and 13-antibody sequencing of CBMCs   human
#> hcabm40k 40,000 Cells From the Human Cell Atlas ICA Bone Marrow Dataset   human
#> ifnb                                  IFNB-Stimulated and Control PBMCs   human
#> panc8                  Eight Pancreas Datasets Across Five Technologies   human
#> pbmc3k                                       3k PBMCs from 10X Genomics   human
#> pbmcsca            Broad Institute PBMC Systematic Comparative Analysis   human
#> stxBrain                        10X Genomics Visium Mouse Brain Dataset   mouse
#>                     system ncells
#> cbmc     CBMC (cord blood)   8617
#> hcabm40k       bone marrow  40000
#> ifnb                  PBMC  13999
#> panc8    Pancreatic Islets  14892
#> pbmc3k                PBMC   2700
#> pbmcsca               PBMC  31021
#> stxBrain             brain  12167
#>                                                                     tech
#> cbmc                                                            CITE-seq
#> hcabm40k                                                          10x v2
#> ifnb                                                              10x v1
#> panc8                   SMARTSeq2, Fluidigm C1, CelSeq, CelSeq2, inDrops
#> pbmc3k                                                            10x v1
#> pbmcsca  10x v2, 10x v3, SMARTSeq2, Seq-Well, inDrops, Drop-seq, CelSeq2
#> stxBrain                                                          visium
#>                                                                              notes
#> cbmc                                                                          <NA>
#> hcabm40k                                                                      <NA>
#> ifnb                                                                          <NA>
#> panc8                                                                         <NA>
#> pbmc3k                                                                        <NA>
#> pbmcsca                                                              HCA benchmark
#> stxBrain One sample split across four datasets as paired anterior/posterior slices
#>          Installed InstalledVersion
#> cbmc          TRUE            3.0.0
#> hcabm40k     FALSE             <NA>
#> ifnb         FALSE             <NA>
#> panc8        FALSE             <NA>
#> pbmc3k        TRUE            3.0.0
#> pbmcsca      FALSE             <NA>
#> stxBrain      TRUE            0.1.1
```

Installation of datasets can be done with `InstallData`; this function will accept either a dataset name (eg. `pbmc3k`) or the corresponding package name (eg. `pbmc3k.SeuratData`). `InstallData` will automatically attach the installed dataset package so one can immediately load and use the dataset.

``` r
InstallData("pbmc3k")
```

Loading a dataset is done using the [`data`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/data.html) function

``` r
data("pbmc3k")
pbmc3k
#> An object of class Seurat 
#> 13714 features across 2700 samples within 1 assay 
#> Active assay: RNA (13714 features, 0 variable features)
```

Dataset documentation and information
-------------------------------------

All datasets provided have help pages built for them. These pages are accessed using the standard [`help`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/help.html) function

``` r
?pbmc3k
```

A full command list for the steps taken to generate each dataset is present in the examples section of these help pages.

Packages will also often have citation information bundled with the package. Citation information can be accessed with `CiteData`

``` r
CiteData("cbmc")
#> 
#> To cite the CBMC dataset, please use:
#> 
#>   Stoeckius et al. Simultaneous epitope and transcriptome measurement
#>   in single cells. Nature Methods (2017)
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     author = {Marlon Stoeckius and Christoph Hafemeister and William Stephenson and Brian Houck-Loomis and Pratip K Chattopadhyay and Harold Swerdlow and Rahul Satija and Peter Smibert},
#>     title = {Simultaneous epitope and transcriptome measurement in single cells},
#>     journal = {Nature Methods},
#>     year = {2017},
#>     doi = {10.1038/nmeth.4380},
#>     url = {https://www.nature.com/articles/nmeth.4380},
#>   }
```

Rationale and Implementation
----------------------------

We created SeuratData in order to distribute datasets for [Seurat](https://satijalab.org/seurat/get_started.html) [vignettes](https://satijalab.org/seurat/frv.html) in as painless and reproducible a way as possible. We also wanted to give users the flexibility to selectively install and load datasets of interest, to minimize disk storage and memory use.

To accomplish this, we opted to distribute datasets through individual R packages. Under the hood, SeuratData uses and extends standard R functions, such as [`install.packages`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/install.packages.html) for dataset installation, [`available.packages`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/available.packages.html) for dataset listing, [`citation`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/citation.html) for dataset citation, and [`data`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/data.html) for dataset loading.

SeuratData therefore serves as a more specific package manager (similar to a metapackage) for R. We provide wrappers around R's package management functions, [extend them to provide relevant metadata](https://github.com/satijalab/seurat-data/#getting-started) about each dataset, and set default settings (for example, the repository where data is stored) to facilitate easy installation.
