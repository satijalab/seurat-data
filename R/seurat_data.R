#' @include zzz.R
#'
NULL

#' Get list of available datasets
#'
#' @return A dataframe with available Seurat datasets. Rownames of the dataframe
#' are the actual package names
#' \describe{
#'   \item{Dataset}{Name of dataset, usable for other functions in SeuratData
#'   (eg. \code{\link{InstallData}})}
#'   \item{Version}{Version of dataset, generally corresponds to version of
#'   Seurat dataset was built with}
#'   \item{Installed}{Is the dataset installed?}
#'   \item{InstalledVersion}{Version of dataset installed, \code{NA}
#'   if not installed}
#' }
#'
#' Other columns are extra metadata, and may change as new datasets
#' are made available
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

#' Get citation information for a dataset
#'
#' @inheritParams utils::citation
#' @param ds Name of dataset to get citation information for
#'
#' @return An object of class \code{citation}, inheriting from class
#' \code{\link[utils]{bibentry}}; see \code{\link[utils]{citation}} for more
#' details
#'
#' @importFrom utils citation
#'
#' @export
#'
#' @seealso \code{\link[utils]{citation}}
#'
#' @examples
#' \dontrun{
#' CiteData('cbmc')
#' }
#'
CiteData <- function(ds, lib.loc = NULL) {
  ds <- ds[1]
  pkg <- NameToPackage(ds = ds)
  return(tryCatch(
    expr = suppressWarnings(expr = citation(
      package = pkg,
      lib.loc = lib.loc,
      auto = FALSE
    )),
    error = function(...) {
      warning(
        "No citation information found for ",
        ds,
        call. = FALSE,
        immediate. = TRUE
      )
      return(invisible(x = NULL))
    }
  ))
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
    installed <- intersect(
      x = pkgs,
      y = NameToPackage(ds = rownames(x = InstalledData()))
    )
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
#' @return A dataframe with installed Seurat datasets. Rownames of the dataframe
#' are the actual package names
#' \describe{
#'   \item{Dataset}{Name of dataset, usable for other functions in SeuratData
#'   (eg. \code{\link{InstallData}})}
#'   \item{Version}{Version of dataset, generally corresponds to version of
#'   Seurat dataset was built with}
#'   \item{Installed}{Is the dataset installed?}
#'   \item{InstalledVersion}{Version of dataset installed}
#' }
#'
#' Other columns are extra metadata, and may change as new datasets
#' are made available
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
#' @inheritParams SeuratDisk::LoadH5Seurat
#' @param ds Optional name of dataset to load
#' @param type How to load the \code{Seurat} object; choose from either
#' 'default' for the default dataset or any dataset listed in the
#' \code{other.datasets} or \code{disk.datasets} sections of the data manifest
#'
#' @return A \code{Seurat} object with the dataset asked for
#'
#' @importFrom utils data
#'
#' @export
#'
#' @seealso \code{\link[utils]{data}}
#'
LoadData <- function(
  ds,
  type = 'default',
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE
) {
  installed <- InstalledData()
  rownames(x = installed) <- NameToPackage(ds = rownames(x = installed))
  if (!NameToPackage(ds = ds) %in% rownames(x = installed)) {
    stop("Cannot find dataset ", ds, call. = FALSE)
  }
  ds <- NameToPackage(ds = ds)
  datasets <- c(
    installed[ds, 'DefaultDS', drop = TRUE],
    trimws(x = unlist(x = strsplit(
      x = installed[ds, 'OtherDS', drop = TRUE],
      split = ','
    ))),
    trimws(x = unlist(x = strsplit(
      x = installed[ds, 'DiskDS', drop = TRUE],
      split = ','
    )))
  )
  type <- match.arg(
    arg = type,
    choices = c('raw', 'default', datasets)
  )
  if (type %in% c('raw', 'default')) {
    type <- gsub(pattern = pkg.key, replacement = '', x = ds)
  } else if (type == 'final') {
    type <- paste0(gsub(pattern = pkg.key, replacement = '', x = ds), '.final')
  }
  if (type %in% data(package = ds)$results[, 'Item', drop = TRUE]) {
    e <- new.env()
    data(list = type, package = ds, envir = e)
    return(e[[type]])
  }
  if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
    stop("Loading a disk-based dataset requires SeuratDisk", call. = FALSE)
  }
  return(SeuratDisk::LoadH5Seurat(
    file = system.file(
      file.path('extdata', 'objects', paste0(type, '.h5Seurat')),
      package = ds,
      mustWork = TRUE
    ),
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    images = images,
    meta.data = meta.data,
    commands = commands,
    misc = misc,
    tools = tools,
    verbose = verbose
  ))
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
