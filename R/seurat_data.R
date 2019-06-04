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
#' @inheritDotParams utils::remove.packages
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
  avail.pkgs$Installed <- vapply(
    X = rownames(x = avail.pkgs),
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE,
  )
  pkg.env$manifest <- avail.pkgs
  invisible(x = NULL)
}
