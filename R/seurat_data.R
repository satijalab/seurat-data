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
  return(pkg.env$manifest)
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
}

#' Get a list of installed datasets
#'
#' @return A vector of installed datasets
#'
#' @export
#'
InstalledData <- function() {
  installed <- vector(mode = 'character')
  for (pkg in rownames(x = pkg.env$manifest)) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      installed <- c(installed, pkg)
    }
  }
  installed <- gsub(pattern = '\\.SeuratData$', replacement = '', x = installed)
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
  pkgs <- pkgs[pkgs %in% pkg.env$manifest$Dataset]
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
  avail.pkgs <- available.packages(
    repos = repo.use,
    type = 'source',
    fields = 'Description'
  )
  avail.pkgs <- as.data.frame(x = avail.pkgs, stringsAsFactors = FALSE)
  avail.pkgs <- subset(
    x = avail.pkgs,
    subset = grepl(pattern = '\\.SeuratData', x = Package)
  )
  avail.pkgs <- apply(
    X = avail.pkgs,
    MARGIN = 1,
    FUN = function(pkg) {
      dataset <- gsub(
        pattern = '\\.SeuratData$',
        replacement = '',
        x = pkg[['Package']]
      )
      # version <- package_version(x = pkg[['Version']])
      desc <- unlist(x = strsplit(x = pkg[['Description']], split = '\n'))
      desc <- sapply(
        X = strsplit(x = desc, split = ':'),
        FUN = function(x) {
          name <- trimws(x = x[[1]])
          val <- trimws(x = unlist(x = strsplit(x = x[[2]], split = ',')))
          val <- paste(val, collapse = ', ')
          names(x = val) <- name
          return(val)
        }
      )
      desc <- sapply(
        X = desc,
        FUN = function(x) {
          x <- tryCatch(
            expr = as.numeric(x = x),
            warning = function(...) {
              return(x)
            }
          )
          if (!is.numeric(x = x) && !is.na(x = as.logical(x = x))) {
            x <- as.logical(x = x)
          }
          return(x)
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
      desc <- c('Dataset' = dataset, 'Version' = pkg[['Version']], desc)
      return(desc)
    }
  )
  avail.pkgs <- lapply(X = avail.pkgs, FUN = as.data.frame, stringsAsFactors = FALSE)
  avail.pkgs <- do.call(what = 'rbind', args = avail.pkgs)
  avail.pkgs$Version <- package_version(x = avail.pkgs$Version)
  pkg.env$manifest <- avail.pkgs
  invisible(x = NULL)
}
