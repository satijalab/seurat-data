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
#' @inheritParams LoadH5Seurat
#' @param ds Optional name of dataset to load
#' @param type How to load the \code{Seurat} object; choose from
#' \describe{
#'   \item{info}{Information about the object and what's stored in it}
#'   \item{raw}{The raw form of the dataset, no other options are evaluated}
#'   \item{processed}{The proccessed data, modular loading avaible by setting other parameters}
#' }
#'
#' @inherit LoadH5Seurat return
#'
#' @importFrom utils data
#'
#' @export
#'
#' @seealso \code{\link[utils]{data}}
#'
LoadData <- function(
  ds,
  type = c('info', 'raw', 'processed'),
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  verbose = TRUE
) {
  .NotYetImplemented()
  installed <- InstalledData()
  if (!NameToPackage(ds = ds) %in% rownames(x = installed)) {
    stop("Cannot find dataset ", ds, call. = FALSE)
  }
  ds <- NameToPackage(ds = ds)
  type <- match.arg(arg = tolower(x = type), choices = c('info', 'raw', 'processed'))
  if (type == 'raw') {
    e <- new.env()
    ds <- gsub(pattern = '\\.SeuratData', replacement = '', x = ds)
    data(list = ds, envir = e)
    return(e[[ds]])
  }
  .NotYetImplemented()
  type <- match.arg(arg = type, choices = c('info', 'processed'))
  return(LoadH5Seurat(
    file = system.file(
      file.path('extdata', 'processed.h5Seurat'),
      package = ds,
      mustWork = TRUE
    ),
    type = ifelse(test = type == 'processed', yes = 'object', no = type),
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    verbose = verbose
  ))
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
