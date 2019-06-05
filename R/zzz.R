#' Install and manage Seurat datasets
#'
#' @docType package
#' @rdname SeuratData-package
#' @name SeuratData-package
#' @keywords internal
#'
"_PACKAGE"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Global variables and environment
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

repo.use <- 'http://satijalab04.nygenome.org'
pkg.env <- new.env()
pkg.env$manifest <- vector(mode = 'character')
pkg.env$attached <- vector(mode = 'character')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Update the available package manifest
#'
#' @return Nothing, updates the manifest in the package environment
#'
#' @importFrom utils available.packages packageVersion
#'
#' @keywords internal
#'
UpdateManifest <- function() {
  avail.pkgs <- available.packages(
    repos = repo.use,
    type = 'source',
    fields = 'Description'
  )
  avail.pkgs <- as.data.frame(x = avail.pkgs, stringsAsFactors = FALSE)
  avail.pkgs <- avail.pkgs[grepl(pattern = '\\.SeuratData$', x = avail.pkgs$Package), , drop = FALSE]
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
  avail.pkgs$InstalledVersion <- sapply(
    X = rownames(x = avail.pkgs),
    FUN = function(x) {
      pkg.version <- tryCatch(
        expr = packageVersion(pkg = x),
        error = function(...) {
          return(NA)
        }
      )
      return(as.character(x = pkg.version))
    }
  )
  avail.pkgs$InstalledVersion <- package_version(
    x = avail.pkgs$InstalledVersion,
    strict = FALSE
  )
  pkg.env$manifest <- avail.pkgs
  invisible(x = NULL)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom utils packageVersion
#'
.onLoad <- function(libname, pkgname) {
  installed <- InstalledData()
  seurat.version <- tryCatch(
    expr = packageVersion(pkg = 'Seurat'),
    error = function(...) {
      return(NA)
    }
  )
  for (pkg in rownames(x = installed)) {
    dataset <- installed[pkg, 'Dataset', drop = TRUE]
    pkg.version <- installed[pkg, 'InstalledVersion', drop = TRUE]
    avail.version <- installed[pkg, 'Version', drop = TRUE]
    if (pkg.version < avail.version) {
      packageStartupMessage(
        "A new version of ",
        dataset,
        " is available, please see UpdateData for help updating your datasets"
      )
    }
    if (!paste0('package', pkg) %in% search()) {
      packageStartupMessage("Attaching ", dataset)
      try(expr = attachNamespace(ns = pkg), silent = TRUE)
      pkg.env$attached <- c(pkg.env$attached, pkg)
    }
    if (!is.na(x = seurat.version) && seurat.version < pkg.version) {
      warning(
        "Dataset ",
        dataset,
        " was built using a newer version of Seurat than currently installed (",
        pkg.version,
        " vs ",
        seurat.version,
        ")",
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
}

.onUnload <- function(libpath) {
  for (pkg in pkg.env$attached) {
    message("Detaching ", pkg.env$manifest[pkg, 'Dataset', drop = TRUE])
    unloadNamespace(ns = pkg)
  }
}
