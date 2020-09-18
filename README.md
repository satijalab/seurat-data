# SeuratData

SeuratData is a mechanism for distributing datasets in the form of [Seurat](https://satijalab.org/seurat) objects using R's internal package and data management systems. It represents an easy way for users to get access to datasets that are used in the Seurat vignettes.

### Installation

Installation of SeuratData can be accomplished through [devtools](https://cran.r-project.org/package=devtools)

```R
devtools::install_github('satijalab/seurat-data')
```

### Getting Started

When loading SeuratData, a list of all available datasets will be displayed (this is similar to other metapackages like [tidyverse](https://cran.r-project.org/package=tidyverse) along with the version of [Seurat](https://satijalab.org/seurat/) used to create each dataset). This message can be suppressed with [`suppressPackageStartupMessages`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/message.html)

```R
> library(SeuratData)
── Installed datasets ───────────────────────────────────────────────────────────── SeuratData v0.1.0 ──
✔ cbmc   3.0.0                                           ✔ panc8  3.0.0
✔ ifnb   3.0.0                                           ✔ pbmc3k 3.0.0

───────────────────────────────────────────────── Key ──────────────────────────────────────────────────
✔ Dataset loaded successfully
```

To see a manifest of all available datasets, use `AvailableData`; this manifest will update as new datasets are uploaded to our data repository.

```R
> AvailableData()
                     Dataset Version                                                        Summary species            system ncells                                                            tech         notes Installed InstalledVersion
cbmc.SeuratData         cbmc   3.0.0                   scRNAseq and 13-antibody sequencing of CBMCs   human CBMC (cord blood)   8617                                                        CITE-seq          <NA>      TRUE            3.0.0
hcabm40k.SeuratData hcabm40k   3.0.0 40,000 Cells From the Human Cell Atlas ICA Bone Marrow Dataset   human       bone marrow  40000                                                          10x v2          <NA>     FALSE            3.0.0
ifnb.SeuratData         ifnb   3.0.0                              IFNB-Stimulated and Control PBMCs   human              PBMC  13999                                                          10x v1          <NA>      TRUE            3.0.0
panc8.SeuratData       panc8   3.0.0               Eight Pancreas Datasets Across Five Technologies   human Pancreatic Islets  14892                SMARTSeq2, Fluidigm C1, CelSeq, CelSeq2, inDrops          <NA>      TRUE            3.0.0
pbmc3k.SeuratData     pbmc3k   3.0.0                                     3k PBMCs from 10X Genomics   human              PBMC   2700                                                          10x v1          <NA>      TRUE            3.0.0
pbmcsca.SeuratData   pbmcsca   3.0.0           Broad Institute PBMC Systematic Comparative Analysis   human              PBMC  31021 10x v2, 10x v3, SMARTSeq2, Seq-Well, inDrops, Drop-seq, CelSeq2 HCA benchmark     FALSE            3.0.0
```

Installation of datasets can be done with `InstallData`; this function will accept either a dataset name (eg. `pbmc3k`) or the corresponding package name (eg. `pbmc3k.SeuratData`). `InstallData` will automatically attach the installed dataset package so one can immediately load and use the dataset.

```R
> InstallData("pbmc3k")
```

Loading a dataset is done using the [`data`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/data.html) function

```R
> data("pbmc3k")
> pbmc3k
An object of class Seurat
13714 features across 2700 samples within 1 assay
Active assay: RNA (13714 features)
```

### Dataset documentation and information

All datasets provided have help pages built for them. These pages are accessed using the standard [`help`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/help.html) function

```R
> ?pbmc3k
> ?ifnb
```

A full command list for the steps taken to generate each dataset is present in the examples section of these help pages.

Packages will also often have citation information bundled with the package. Citation information can be accessed by passing the **package name**, not the dataset name, to the [`citation`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/citation.html) function

```R
> citation('cbmc.SeuratData')

To cite the CBMC dataset, please use:

  Stoeckius et al. Simultaneous epitope and transcriptome measurement in
  single cells. Nature Methods (2017)

A BibTeX entry for LaTeX users is

  @Article{,
    author = {Marlon Stoeckius and Christoph Hafemeister and William Stephenson and Brian Houck-Loomis and Pratip K Chattopadhyay and Harold Swerdlow and Rahul Satija and Peter Smibert},
    title = {Simultaneous epitope and transcriptome measurement in single cells},
    journal = {Nature Methods},
    year = {2017},
    doi = {10.1038/nmeth.4380},
    url = {https://www.nature.com/articles/nmeth.4380},
  }
```

### Rationale and Implementation

We created SeuratData in order to distribute datasets for [Seurat](https://satijalab.org/seurat/get_started.html) [vignettes](https://satijalab.org/seurat/frv.html) in as painless and reproducible a way as possible. We also wanted to give users the flexibility to selectively install and load datasets of interest, to minimize disk storage and memory use.

To accomplish this, we opted to distribute datasets through individual R packages. Under the hood, SeuratData uses and extends standard R functions, such as [`install.packages`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/install.packages.html) for dataset installation, [`available.packages`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/available.packages.html) for dataset listing, and [`data`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/data.html) for dataset loading. 

SeuratData therefore serves as a more specific package manager (similar to a metapackage) for R. We provide wrappers around R's package management functions, [extend them to provide relevant metadata](https://github.com/satijalab/seurat-data/#getting-started) about each dataset, and set default settings (for example, the repository where data is stored) to facilitate easy installation.
