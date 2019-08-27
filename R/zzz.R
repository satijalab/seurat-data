#' Install and manage Seurat datasets
#'
#' @section Package options:
#'
#' SeuratData uses the following [options()] to configure behaviour:
#'
#' \itemize{
#'   \item `SeuratData.repo.use`: Set the location where the SeuratData datasets
#'   are stored. Users generally should not modify.
#' }
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

default.options <- list(
  SeuratData.repo.use = 'https://seurat.nygenome.org'
)

pkg.env <- new.env()
pkg.env$manifest <- vector(mode = 'character')
pkg.env$attached <- vector(mode = 'character')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Set a default value if an object is null
#'
#' @param lhs An object to set if it's null
#' @param rhs The value to provide if x is null
#'
#' @return rhs if lhs is null, else lhs
#'
#' @name set-if-null
#'
#' @author Hadley Wickham
#' @references \url{https://adv-r.hadley.nz/functions.html#missing-arguments}
#'
#' @examples
#' \dontrun{
#' 4 %||% 5
#' NULL %|| 5
#' }
#'
#' @keywords internal
#'
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

#' Attach datasets
#'
#' @param pkgname Name of package
#'
#' @return Invisible \code{NULL}
#'
#' @importFrom stats na.omit
#' @importFrom cli rule symbol
#' @importFrom utils packageVersion
#' @importFrom crayon bold red green yellow blue col_align col_nchar
#'
#' @keywords internal
#'
AttachData <- function(pkgname = 'SeuratData') {
  installed <- InstalledData()
  installed <- installed[!paste0('package:', rownames(x = installed)) %in% search(), , drop = FALSE]
  installed <- installed[grep(pattern = 'NA', x = rownames(x = installed), invert = TRUE), ]
  if (nrow(x = installed) < 1) {
    return(invisible(x = NULL))
  }
  seurat.version <- tryCatch(
    expr = packageVersion(pkg = 'Seurat'),
    error = function(...) {
      return(NA)
    }
  )
  header <- rule(
    left = bold('Installed datasets'),
    right = paste0(pkgname, ' v', packageVersion(pkg = pkgname))
  )
  symbols <- if (is.na(x = seurat.version)) {
    red(symbol$fancy_question_mark)
  } else {
    c(green(symbol$tick), yellow(symbol$pointer))[(installed$InstalledVersion > seurat.version) + 1]
  }
  packageStartupMessage(header)
  pkgs <- paste(
    symbols,
    blue(format(x = installed$Dataset)),
    col_align(text = installed$InstalledVersion, width = max(col_nchar(x = installed$InstalledVersion)))
  )
  if (length(x = pkgs) %% 2 == 1) {
    pkgs <- c(pkgs, '')
  }
  col.1 <- seq_len(length.out = length(x = pkgs) / 2)
  space.start <- floor(min(
    col_nchar(x = header) - max(col_nchar(x = pkgs[-col.1])),
    col_nchar(x = header) * (1 / 2)
  ))
  space <- MakeSpace(n = max(space.start - max(col_nchar(x = pkgs[col.1])), 1))
  info <- paste0(pkgs[col.1], space, pkgs[-col.1])
  packageStartupMessage(paste(info, collapse = '\n'), '\n')
  packageStartupMessage(rule(center = 'Key'))
  packageStartupMessage(
    paste(
      c(green(symbol$tick), yellow(symbol$pointer), red(symbol$fancy_question_mark)),
      c(
        'Dataset loaded successfully',
        'Dataset built with a newer version of Seurat than installed',
        'Unknown version of Seurat installed'
      ),
      collapse = '\n'
    ),
    '\n'
  )
  suppressPackageStartupMessages(expr = lapply(
    X = rownames(x = installed),
    FUN = attachNamespace
  ))
  pkg.env$attached <- rownames(x = installed)
  invisible(x = NULL)
}

#' Enumerate a list or vector
#'
#' @param x A list or a vector
#'
#' @return A list of length \code{x} with the following named values:
#' \describe{
#'   \item{\code{name}}{The name of \code{x} at a given index}
#'   \item{\code{value}}{The value of \code{x} at the corresponding index}
#' }
#'
#' @note For any given index \code{i} in \code{x}, all attempts to use the name
#' of the value of \code{x} at \code{i} will be made. If there is no name
#' (eg. \code{nchar(x = names(x = x)[i]) == 0}), the index \code{i} will be used
#' in its stead
#'
#' @keywords internal
#'
Enumerate <- function(x) {
  indices <- seq_along(along.with = x)
  keys <- names(x = x) %||% as.character(x = indices)
  keys[nchar(x = keys) == 0] <- indices[nchar(x = keys) == 0]
  vals <- lapply(
    X = indices,
    FUN = function(i) {
      return(c('name' = keys[i], 'value' = unname(obj = x[i])))
    }
  )
  return(vals)
}

#' Check to see if a matrix is empty
#'
#' Deterime if a matrix is empty or not
#'
#' @param x A matrix
#'
#' @return \code{TRUE} if the matrix is empty otherwise \code{FALSE}
#'
#' @details A matrix is considered empty if it satisfies one of the following
#' conditions:
#' \itemize{
#'   \item The dimensions of the matrix are 0-by-0 (\code{all(dim(x) == 0)})
#'   \item The dimensions of the matrix are 1-by-1 (\code{all(dim(x) == 1)}) and
#'   the sole value is a \code{NA}
#' }
#' These two situations correspond to matrices generated with either
#' \code{new('matrix')} or \code{matrix()}
#'
#' @examples
#' \dontrun{
#' IsMatrixEmpty(new('matrix'))
#' IsMatrixEmpty(matrix())
#' IsMatrixEmpty(matrix(1:9, nrow = 3))
#' }
#'
#' @keywords internal
#'
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

#' Make a space
#'
#' @param n Length space should be
#'
#' @return A space (' ') of length \code{n}
#'
#' @examples
#' \dontrun{
#' MakeSpace(10)
#' }
#'
#' @keywords internal
#'
MakeSpace <- function(n) {
  return(paste(rep_len(x = ' ', length.out = n), collapse = ''))
}

#' Find dataset packages based on name
#'
#' @param ds Names of datasets
#'
#' @return A vector of package names based on dataset names
#'
#' @examples
#' \dontrun{
#' NameToPackage('cbmc')
#' NameToPackage('pbmc3k.SeuratData')
#' NameToPackage('notadataset')
#' }
#'
#' @keywords internal
#'
NameToPackage <- function(ds) {
  avail.pkgs <- AvailableData()
  ds.use <- ds %in% c(rownames(x = avail.pkgs), avail.pkgs$Dataset)
  if (sum(x = ds.use) != length(x = ds)) {
    warning(
      "The following datasets are not available: ",
      paste(ds[!ds.use], collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  ds <- ds[ds.use]
  ds.replace <- which(x = ds %in% avail.pkgs$Dataset)
  pkg.replace <- which(x = avail.pkgs$Dataset %in% ds[ds.replace])
  ds[ds.replace] <- rownames(x = avail.pkgs)[pkg.replace]
  if (length(x = ds) < 1) {
    stop("No SeuratData datasets provided", call. = FALSE)
  }
  return(ds)
}

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
    repos = getOption(x = "SeuratData.repo.use"),
    type = 'source',
    fields = c('Description', 'Title'),
    ignore_repo_cache = TRUE
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
      desc <- c(
        'Dataset' = dataset,
        'Version' = pkg[['Version']],
        'Summary' = pkg[['Title']],
        desc
      )
      return(desc)
    }
  )
  manifest.names <- unique(x = unlist(
    x = lapply(X = avail.pkgs, FUN = names),
    use.names = FALSE
  ))
  for (pkg in names(x = avail.pkgs)) {
    for (col in manifest.names) {
      avail.pkgs[[pkg]][[col]] <- avail.pkgs[[pkg]][[col]] %||% NA
    }
  }
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
  # # TODO: remove these when we allow loading of processed datasets
  # ds.index <- which(x = colnames(x = avail.pkgs) %in% c('default.dataset', 'other.datasets'))
  # avail.pkgs <- avail.pkgs[, -ds.index]
  pkg.env$manifest <- avail.pkgs
  invisible(x = NULL)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onAttach <- function(libname, pkgname) {
  AttachData(pkgname = pkgname)
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(x = default.options) %in% names(x = op))
  if (any(toset)) {
    options(default.options[toset])
  }
  invisible()
}
