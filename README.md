# SeuratData

SeuratData is a mechanism for distributing datasets in the form of [Seurat](https://satijalab.org/seurat) objects using R's internal package and data management systems.

### Rationale and Implementation

Dataset distribution is something that should be done in a standard and, hopefully, painless manner. The goal for SeuratData is to distribute datasets for [Seurat](https://satijalab.org/seurat/get_started.html) [vignettes](https://satijalab.org/seurat/frv.html) in a manner that doesn't depend downloading files from various locations across the internet, reading them into R, and keeping track of them on the filesystem for the life of the vignettes. We also wanted to selectively distribute datasets, so not all datasets are downloaded and brought into memory at once, in a manner that was easy to update and as R-native as possible.

To accomplish this, we opted to distribute datasets through individual R packages. Under the hood, SeuratData uses standard R functions, such as [`install.packages`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/install.packages.html) for dataset installation, [`available.packages`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/available.packages.html) for dataset listing, and [`data`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/data.html) for dataset loading. However, these functions are fairly generalized and provide little useful metadata for datasets.

SeuratData serves as a more specific package manager and quasi-metapacakge for R, designed around making datasets easy to install and get information for. We provide wrappers around R's package management functions with sensible defaults for dataset packages, and [extend them to provide relavent metadata](https://github.com/satijalab/seurat-data/#getting-started) about said dataset packages. These defaults include setting the package repository (dataset packages are generally not hosted on current R package repositories, and submitting them there would involve long waits before new data is available; to see the URL for the dataset repository, use `getOption("SeuratData.repo.use")`) and building the datasets from source (datasets have no compiled code, and therefore don't need binary installations). Finally, SeuratData will automatically attach dataset packages upon loading, allowing users to load the datasets into memory without needing to load the dataset packages individually.

### Installation

Installation of SeuratData can be accomplished through [devtools](https://cran.r-project.org/package=devtools)

```R
devtools::install_github('satijalab/seurat-data')
```

### Getting Started

When loading SeuratData, any and all datasets installed through SeuratData will automatically be attached as well. Similar to other metapackages like [tidyverse](https://cran.r-project.org/package=tidyverse), a list of attached datasets will be displayed. Dataset versions correspond to the version of [Seurat](https://satijalab.org/seurat/) they were built under. This message can be suppressed with [`suppressPackageStartupMessages`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/message.html)

```R
> library(SeuratData)
── Attaching datasets ──────────────────────────────────────────────────────────── SeuratData v0.1.0 ──
✔ cbmc   3.0.0                                           ✔ panc8  3.0.0
✔ ifnb   3.0.0                                           ✔ pbmc3k 3.0.0

───────────────────────────────────────────────── Key ──────────────────────────────────────────────────
✔ Dataset loaded succesfully
❯ Dataset built with a newer version of Seurat than installed
❓ Unknown version of Seurat installed

```

To see a manifest of available datasets, use `AvailableData`; this manifest will update as new datasets are uploaded to our data repository.

```R
> AvailableData()
                     Dataset Version                                                        Summary species            system ncells                                                            tech default.dataset other.datasets         notes Installed InstalledVersion
cbmc.SeuratData         cbmc   3.0.0                   scRNAseq and 13-antibody sequencing of CBMCs   human CBMC (cord blood)   8617                                                        CITE-seq             raw           <NA>          <NA>      TRUE            3.0.0
hcabm40k.SeuratData hcabm40k   3.0.0 40,000 Cells From the Human Cell Atlas ICA Bone Marrow Dataset   human       bone marrow  40000                                                          10x v2             raw           <NA>          <NA>     FALSE            3.0.0
ifnb.SeuratData         ifnb   3.0.0                              IFNB-Stimulated and Control PBMCs   human              PBMC  13999                                                          10x v1             raw           <NA>          <NA>      TRUE            3.0.0
panc8.SeuratData       panc8   3.0.0               Eight Pancreas Datasets Across Five Technologies   human Pancreatic Islets  14892                SMARTSeq2, Fluidigm C1, CelSeq, CelSeq2, inDrops             raw           <NA>          <NA>      TRUE            3.0.0
pbmc3k.SeuratData     pbmc3k   3.0.0                                     3k PBMCs from 10X Genomics   human              PBMC   2700                                                          10x v1             raw          final          <NA>      TRUE            3.0.0
pbmcsca.SeuratData   pbmcsca   3.0.0           Broad Institute PBMC Systematic Comparative Analysis   human              PBMC  31021 10x v2, 10x v3, SMARTSeq2, Seq-Well, inDrops, Drop-seq, CelSeq2             raw           <NA> HCA benchmark     FALSE            3.0.0
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

Using just the dataset name refers to the "default dataset"; this is generally a raw version of the dataset with no processing of the data. Certain packages have multiple forms of the dataset, such as a processed or "final" version. Loading these alternate versions is done by appending the modifier to the dataset name with a period

```R
> data("pbmc3k.final")
> pbmc3k.final
An object of class Seurat
13714 features across 2638 samples within 1 assay
Active assay: RNA (13714 features)
 2 dimensional reductions calculated: pca, umap
```

### Dataset documentation and information

All datasets provided have help pages built for them. These pages are accessed using the standard [`help`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/help.html) function

```R
> ?pbmc3k
> ?pbmc3k.final
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
