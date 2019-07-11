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
#' @seealso \code{\link{InstallData}} \code{\link{InstalledData}} \code{\link{RemoveData}} \code{\link{UpdateData}}
#'
AvailableData <- function() {
  UpdateManifest()
  return(pkg.env$manifest)
}

#' Install a dataset
#'
#' @param ds List of datasets to install
#' @param ... Extra parameters passed to \code{\link[utils]{install.packages}}
#'
#' @return Invisible \code{NULL}
#'
#' @importFrom utils install.packages
#'
#' @export
#'
#' @seealso \code{\link{AvailableData}} \code{\link{InstalledData}} \code{\link{RemoveData}} \code{\link{UpdateData}}
#'
InstallData <- function(ds, ...) {
  UpdateManifest()
  pkgs <- NameToPackage(ds = ds)
  pkgs2 <- paste0('package:', pkgs)
  for (p in pkgs2[pkgs2 %in% search()]) {
    detach(name = p, unload = TRUE, character.only = TRUE)
  }
  install.packages(pkgs = pkgs, repos = getOption(x = "SeuratData.repo.use"), type = 'source', ...)
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
#' @seealso \code{\link{AvailableData}} \code{\link{InstallData}} \code{\link{RemoveData}} \code{\link{UpdateData}}
#'
InstalledData <- function() {
  dat <- AvailableData()
  return(dat[which(x = dat$Installed, ), , drop = FALSE])
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
#' @seealso \code{\link{AvailableData}} \code{\link{InstallData}} \code{\link{InstalledData}} \code{\link{UpdateData}}
#'
RemoveData <- function(ds, lib) {
  UpdateManifest()
  pkgs <- NameToPackage(ds = ds)
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
#' @seealso \code{\link{AvailableData}} \code{\link{InstallData}} \code{\link{InstalledData}} \code{\link{RemoveData}}
#'
UpdateData <- function(ask = TRUE, lib.loc = NULL) {
  update.packages(lib.loc = lib.loc, repos = getOption(x = "SeuratData.repo.use"), ask = ask, type = 'source')
  UpdateManifest()
  invisible(x = NULL)
}
