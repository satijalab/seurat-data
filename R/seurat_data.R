#' @include zzz.R
#'
NULL

#' Get list of available datasets
#'
#' @return A dataframe with available Seurat datasets. Rownames of the dataframe are the actual package names
#' \describe{
#'   \item{Dataset}{Name of dataset, usable for other functions in SeuratData (eg. \code{\link{InstallData}})}
#'   \item{Version}{Version of dataset, generally corresponds to version of Seurat dataset was built with}
#'   \item{Installed}{Is the dataset installed?}
#'   \item{InstalledVersion}{Version of dataset installed, \code{NA} if not installed}
#' }
#'
#' Other columns are extra metadata, and may change as new datasets are made available
#'
#' @export
#'
#' @seealso \code{\link{InstallData}} \code{\link{InstalledData}}
#' \code{\link{RemoveData}} \code{\link{UpdateData}}
#'
#' @examples
#' AvailableData()
#'
AvailableData <- function() {
  UpdateManifest()
  return(pkg.env$manifest)
}

#' Install a dataset
#'
#' @param ds List of datasets to install
#' @param force.reinstall Force reinstall already installed packages
#' @param ... Extra parameters passed to \code{\link[utils]{install.packages}}
#'
#' @return Invisible \code{NULL}
#'
#' @importFrom utils install.packages
#'
#' @export
#'
#' @seealso \code{\link{AvailableData}} \code{\link{InstalledData}}
#' \code{\link{RemoveData}} \code{\link{UpdateData}}
#'
#' @examples
#' \dontrun{
#' InstallData('pbmc3k')
#' }
#'
InstallData <- function(ds, force.reinstall = FALSE, ...) {
  UpdateManifest()
  if (pkg.env$source != 'remote') {
    stop(
      "No access to remote SeuratData repository, unable to install new datasets",
      call. = FALSE
    )
  }
  pkgs <- NameToPackage(ds = ds)
  if (!force.reinstall) {
    installed <- intersect(x = pkgs, y = rownames(x = InstalledData()))
    if (length(x = installed) > 0) {
      warning(
        "The following packages are already installed and will not be reinstalled: ",
        paste(
          gsub(pattern = pkg.key, replacement = '', x = installed),
          collapse = ', '
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      pkgs <- setdiff(x = pkgs, y = installed)
      if (!as.logical(x = length(x = pkgs))) {
        UpdateManifest()
        return(invisible(x = NULL))
      }
    }
  }
  pkgs2 <- paste0('package:', pkgs)
  for (p in pkgs2[pkgs2 %in% search()]) {
    detach(name = p, unload = TRUE, character.only = TRUE)
  }
  install.packages(
    pkgs = pkgs,
    repos = getOption(x = "SeuratData.repo.use"),
    type = 'source',
    ...
  )
  for (pkg in pkgs) {
    attachNamespace(ns = pkg)
    pkg.env$attached <- c(pkg.env$attached, pkg)
  }
  UpdateManifest()
  invisible(x = NULL)
}

#' Get a list of installed datasets
#'
#' @return A dataframe with installed Seurat datasets. Rownames of the dataframe are the actual package names
#' \describe{
#'   \item{Dataset}{Name of dataset, usable for other functions in SeuratData (eg. \code{\link{InstallData}})}
#'   \item{Version}{Version of dataset, generally corresponds to version of Seurat dataset was built with}
#'   \item{Installed}{Is the dataset installed?}
#'   \item{InstalledVersion}{Version of dataset installed}
#' }
#'
#' Other columns are extra metadata, and may change as new datasets are made available
#'
#' @export
#'
#' @seealso \code{\link{AvailableData}} \code{\link{InstallData}}
#' \code{\link{RemoveData}} \code{\link{UpdateData}}
#'
#' @examples
#' InstalledData()
#'
InstalledData <- function() {
  dat <- AvailableData()
  return(dat[which(x = dat$Installed, ), , drop = FALSE])
}

#' Modularly load a dataset
#'
# @inheritParams LoadH5Seurat
#' @param ds Optional name of dataset to load
#' @param type How to load the \code{Seurat} object; choose from either
#' 'default' for the default dataset or any dataset listed in the
#' \code{other.datasets} section of the data manifest
# \describe{
#   \item{info}{Information about the object and what's stored in it}
#   \item{raw}{The raw form of the dataset, no other options are evaluated}
#   \item{azimuth}{Load the dataset as an azimuth reference}
#   \item{processed}{The proccessed data, modular loading avaible by setting other parameters}
# }
#'
# @inherit LoadH5Seurat return
#' @return A \code{Seurat} object with the dataset asked for
#'
#' @importFrom utils data
#' @importFrom SeuratObject Assays UpdateSeuratObject
#' @importFrom methods as
#'
#' @export
#'
#' @seealso \code{\link[utils]{data}}
#'
LoadData <- function(
  ds,
  type = 'default'
  # assays = NULL,
  # reductions = NULL,
  # graphs = NULL,
  # verbose = TRUE
) {
  installed <- InstalledData()
  if (!NameToPackage(ds = ds) %in% rownames(x = installed)) {
    stop("Cannot find dataset ", ds, call. = FALSE)
  }
  ds <- NameToPackage(ds = ds)
  datasets <- c(
    installed[ds, 'default.dataset', drop = TRUE],
    trimws(x = unlist(x = strsplit(
      x = installed[ds, 'other.datasets', drop = TRUE], split = ','
    )))
  )
  type <- match.arg(
    arg = tolower(x = type),
    choices = c('raw', 'default', 'azimuth', datasets)
  )
  if (type == 'azimuth') {
    return(Azimuth::LoadReference(file.path(system.file(package=ds), type)))
  } else if (type %in% c('raw', 'default')) {
    type <- gsub(pattern = pkg.key, replacement = '', x = ds)
  } else if (type == 'final') {
    type <- paste0(gsub(pattern = pkg.key, replacement = '', x = ds), '.final')
  }
  if (type %in% data(package = ds)$results[, 'Item', drop = TRUE]) {
    e <- new.env()
    data(list = type, package = ds, envir = e)
    # ds <- gsub(pattern = '\\.SeuratData', replacement = '', x = ds)
    # data(list = ds, envir = e)
    e[[type]] <- UpdateSeuratObject(e[[type]])
    assay_option <- getOption(
      x = 'Seurat.object.assay.version',
      default =  'v5'
    )
    if (assay_option == 'v5') {
      update.assays <- intersect(Assays(e[[type]]), c('RNA', 'ADT'))
      for (i in seq_along(along.with = update.assays)) {
        suppressPackageStartupMessages(
          expr = e[[type]][[update.assays[i]]] <- as(
          object = e[[type]][[update.assays[i]]],
          Class = 'Assay5'
          )
        )
      }
    }
    return(e[[type]])
  }
  stop(
    "Could not find dataset '",
    type,
    "', please check manifest and try again",
    call. = FALSE
  )
  .NotYetImplemented()
  # type <- match.arg(arg = type, choices = c('info', 'processed'))
  # return(LoadH5Seurat(
  #   file = system.file(
  #     file.path('extdata', 'processed.h5Seurat'),
  #     package = ds,
  #     mustWork = TRUE
  #   ),
  #   type = ifelse(test = type == 'processed', yes = 'object', no = type),
  #   assays = assays,
  #   reductions = reductions,
  #   graphs = graphs,
  #   verbose = verbose
  # ))
}

#' Remove a dataset
#'
#' @inheritParams utils::remove.packages
#' @param ds List of datasets to remove
#'
#' @importFrom utils remove.packages
#'
#' @export
#'
#' @seealso \code{\link{AvailableData}} \code{\link{InstallData}}
#' \code{\link{InstalledData}} \code{\link{UpdateData}}
#'
#' @examples
#' \dontrun{
#' RemoveData('pbmc3k')
#' }
#'
RemoveData <- function(ds, lib) {
  UpdateManifest()
  pkgs <- NameToPackage(ds = ds)
  pkgs2 <- paste0('package:', pkgs)
  for (p in pkgs2[pkgs2 %in% search()]) {
    detach(name = p, unload = TRUE, character.only = TRUE)
  }
  remove.packages(pkgs = pkgs, lib = lib)
  UpdateManifest()
}

#' Update datasets
#'
#' @inheritParams utils::update.packages
#'
#' @return Invisible \code{NULL}
#'
#' @importFrom utils update.packages
#'
#' @export
#'
#' @seealso \code{\link{AvailableData}} \code{\link{InstallData}}
#' \code{\link{InstalledData}} \code{\link{RemoveData}}
#'
#' @examples
#' \dontrun{
#' UpdateData(ask = FALSE)
#' }
#'
UpdateData <- function(ask = TRUE, lib.loc = NULL) {
  UpdateManifest()
  if (pkg.env$source != 'remote') {
    stop(
      "No access to remote SeuratData repository, unable to update datasets",
      call. = FALSE
    )
  }
  update.packages(lib.loc = lib.loc, repos = getOption(x = "SeuratData.repo.use"), ask = ask, type = 'source')
  UpdateManifest()
  invisible(x = NULL)
}

#' Load the reference RDS files
#'
#' Read in a reference \code{\link[Seurat]{Seurat}} object and annoy index. This
#' function can read from a file path. In order to read properly,
#' there must be the following files:
#' \itemize{
#'  \item \dQuote{ref.Rds} for the downsampled reference \code{Seurat}
#'  object (for mapping)
#'  \item \dQuote{idx.annoy} for the nearest-neighbor index object
#' }
#'
#' @param path Pathto the two RDS files
#'
#' @return A list with two entries:
#' \describe{
#'  \item{\code{map}}{
#'   The downsampled reference \code{\link[Seurat]{Seurat}}
#'   object (for mapping)
#'  }
#'  \item{\code{plot}}{The reference \code{Seurat} object (for plotting)}
#' }
#'
#' @importFrom SeuratObject Idents<- Tool Cells Misc AddMetaData CreateSeuratObject
#' @importFrom Seurat LoadAnnoyIndex
#' @importFrom utils download.file
#' @importFrom Matrix sparseMatrix
#' @importFrom methods slot
#'
LoadReference <- function(path) {
  ref.names <- list(
    map = 'ref.Rds',
    ann = 'idx.annoy'
  )
  mapref <- file.path(path, ref.names$map)
  annref <- file.path(path, ref.names$ann)
  exists <- file.exists(c(mapref, annref))
  if (!all(exists)) {
    stop(
      "Missing the following files from the directory provided: ",
      unlist(x = ref.names)[!exists]
    )
  }
  # Load the map reference
  map <- readRDS(file = mapref)
  # Load the annoy index into the Neighbor object in the neighbors slot
  map[["refdr.annoy.neighbors"]] <- LoadAnnoyIndex(
    object = map[["refdr.annoy.neighbors"]],
    file = annref
  )
  # Create plotref
  ad <- Tool(object = map, slot = "AzimuthReference")
  # plotref.dr <- GetPlotRef(object = ad)
  print(ad)
  plotref.dr <- slot(object = ad, name = "plotref")
  cm <- sparseMatrix(
    i = 1, j = 1, x = 0, dims = c(1, nrow(x = plotref.dr)),
    dimnames = list("placeholder", Cells(x = plotref.dr))
  )
  plot <- CreateSeuratObject(
    counts = cm
  )
  plot[["refUMAP"]] <- plotref.dr
  plot <- AddMetaData(object = plot, metadata = Misc(object = plotref.dr, slot = "plot.metadata"))
  gc(verbose = FALSE)
  return(list(
    map = map,
    plot = plot
  ))
}
