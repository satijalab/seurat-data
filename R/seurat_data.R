#' @include zzz.R
#'
NULL

#' Get list of available datasets
#'
#' @return A vector of dataset names
#'
#' @export
#'
AvailableData <- function() {
  UpdateManifest()
  return(manifest)
}

#' Install a dataset
#'
#' @param pkgs List of datasets to install
#' @param ... Extra parameters passed to \code{\link[utils]{install.packages}}
#'
#' @importFrom utils install.packages
#'
#' @export
#'
InstallData <- function(pkgs, ...) {
  UpdateManifest()
  pkgs <- pkgs[pkgs %in% manifest]
  if (length(x = pkgs) < 1) {
    stop("No datasets provided", call. = FALSE)
  }
  install.packages(
    pkgs = pkgs,
    repos = 'http://satijalab04.nygenome.org',
    type = 'source',
    ...
  )
  for (pkg in pkgs) {
    attachNamespace(ns = pkg)
    attached <<- c(attached, pkg)
  }
}

#' Get a list of installed datasets
#'
#' @return A vector of installed datasets
#'
#' @export
#'
InstalledData <- function() {
  installed <- vector(mode = 'character')
  for (pkg in manifest) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      installed <- c(installed, pkg)
    }
  }
  return(installed)
}

#' Remove a dataset
#'
#' @inheritDotParams utils::remove.packages
#' @param pkgs List of dataset to remove
#'
#' @importFrom utils remove.packages
#'
RemoveData <- function(pkgs, lib) {
  UpdateManifest()
  pkgs <- pkgs[pkgs %in% manifest]
  if (length(x = pkgs) < 1) {
    stop("No datasets provided", call. = FALSE)
  }
  remove.packages(pkgs = pkgs, lib = lib)
}

#' Update the available package manifest
#'
#' @importFrom utils available.packages
#'
UpdateManifest <- function() {
  repo.use <- 'http://satijalab04.nygenome.org'
  avail.pkgs <- available.packages(
    contriburl = paste(repo.use, 'src/contrib', sep = '/'),
    type = 'source'
  )
  manifest <<- rownames(x = avail.pkgs)
  invisible(x = NULL)
}