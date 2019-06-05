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
#' @seealso \code{\link{InstallData}} \code{\link{InstalledData}} \code{\link{RemoveData}}
#'
AvailableData <- function() {
  UpdateManifest()
  return(pkg.env$manifest)
}

#' Install a dataset
#'
#' @param pkgs List of datasets to install
#' @param ... Extra parameters passed to \code{\link[utils]{install.packages}}
#'
#' @return Invisible \code{NULL}
#'
#' @importFrom utils install.packages
#'
#' @export
#'
#' @seealso \code{\link{AvailableData}} \code{\link{InstalledData}} \code{\link{RemoveData}}
#'
InstallData <- function(pkgs, ...) {
  UpdateManifest()
  pkgs <- pkgs[pkgs %in% pkg.env$manifest$Dataset]
  if (length(x = pkgs) < 1) {
    stop("No datasets provided", call. = FALSE)
  }
  pkgs <- paste0(pkgs, '.SeuratData')
  install.packages(
    pkgs = pkgs,
    repos = repo.use,
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
#' @seealso \code{\link{AvailableData}} \code{\link{InstallData}} \code{\link{RemoveData}}
#'
InstalledData <- function() {
  dat <- AvailableData()
  return(dat[which(x = dat$Installed, ), , drop = FALSE])
}

#' Remove a dataset
#'
#' @inheritParams utils::remove.packages
#' @param pkgs List of datasets to remove
#'
#' @importFrom utils remove.packages
#'
#' @export
#'
#' @seealso \code{\link{AvailableData}} \code{\link{InstallData}} \code{\link{InstalledData}}
#'
RemoveData <- function(pkgs, lib) {
  UpdateManifest()
  pkgs <- pkgs[pkgs %in% pkg.env$manifest$Dataset]
  if (length(x = pkgs) < 1) {
    stop("No datasets provided", call. = FALSE)
  }
  remove.packages(pkgs = pkgs, lib = lib)
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
UpdateData <- function(ask = TRUE, lib.loc = NULL) {
  update.packages(lib.loc = lib.loc, repos = repo.use, ask = ask, type = 'source')
  UpdateManifest()
  invisible(x = NULL)
}
