#' Install and manage Seurat datasets
#'
#' @section Package options:
#'
#' SeuratData uses the following options to control behaviour, users can configure
#' these with \code{\link[base]{options}}:
#'
#' \itemize{
#'   \item `SeuratData.repo.use`: Set the location where the SeuratData datasets
#'   are stored. Users generally should not modify.
#'   \item `SeuratData.manifest.cache`: Cache the data manifest whenever we talk
#'   to the data repository; note, setting to \code{FALSE} will simply prevent
#'   SeuratData from caching the manifest, not from reading a previously cached
#'   manifest
#'   \item `SeuratData.roaming`: For Windows users, use a roaming profile directory
#'   for domain users. See \url{https://en.wikipedia.org/wiki/Roaming_user_profile}
#'   for a brief overview and Microsoft's documentation for greater detail
#' }
#'
#' @docType package
#' @rdname SeuratData-package
#' @name SeuratData-package
#' @keywords internal
#'
"_PACKAGE"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Global variables, default options, and package environment
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

default.options <- list(
  SeuratData.repo.use = 'http://seurat.nygenome.org/',
  SeuratData.manifest.cache = TRUE,
  SeuratData.roaming = FALSE
)

pkg.key <- '\\.SeuratData$'

pkg.env <- new.env()
pkg.env$manifest <- vector(mode = 'list')
pkg.env$source <- vector(mode = 'character')
pkg.env$update.call <- vector(mode = 'character')
pkg.env$attached <- vector(mode = 'character')
pkg.env$extdata.warn <- FALSE

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

#' Get an application data directory
#'
#' @inheritParams AttachData
#' @param author Author name for application
#'
#' @return A character vector with path to the application data dir
#'
#' @importFrom rappdirs user_data_dir
#'
#' @seealso \code{\link[rappdirs]{user_data_dir}}
#'
#' @keywords internal
#'
AppData <- function(pkgname = 'SeuratData', author = pkgname) {
  return(user_data_dir(
    appname = pkgname,
    appauthor = author,
    version = file.path('%p-platform', '%v'),
    roaming = getOption(x = 'SeuratData.roaming', default = FALSE)
  ))
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
  space <- paste(
    rep_len(
      x = ' ',
      length.out = max(space.start - max(col_nchar(x = pkgs[col.1])), 1)
    ),
    collapse = ''
  )
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
  # Set some defaults
  pkg.env$source <- 'remote'
  cache.manifest <- file.path(
    AppData(pkgname = 'SeuratData', author = 'Satija Lab'),
    'manifest.Rds'
  )
  # Attempt to get the manifest from the remote server
  avail.pkgs <- tryCatch(
    expr = available.packages(
      repos = getOption(x = "SeuratData.repo.use"),
      type = 'source',
      fields = c('Description', 'Title'),
      ignore_repo_cache = TRUE
    ),
    warning = function(...) {
      pkg.env$source <- ifelse(
        test = file.exists(cache.manifest),
        yes = 'appdir',
        no = 'extdata'
      )
    }
  )
  # Process the manifest
  if (pkg.env$source == 'remote') {
    # Lots of stuff to get the manifest modified from the available.packages format
    # into something usable by SeuratData
    pkg.env$extdata.warn <- FALSE
    avail.pkgs <- as.data.frame(x = avail.pkgs, stringsAsFactors = FALSE)
    # Ensure we only use datasets tagged with .SeuratData
    avail.pkgs <- avail.pkgs[grepl(pattern = pkg.key, x = avail.pkgs$Package), , drop = FALSE]
    # Filter down to dataset name, short summary from package title, and
    # metadata contained in package description
    avail.pkgs <- apply(
      X = avail.pkgs,
      MARGIN = 1,
      FUN = function(pkg) {
        # Get dataset name
        dataset <- gsub(
          pattern = pkg.key,
          replacement = '',
          x = pkg[['Package']]
        )
        # Process the description metadata
        desc <- pkg[['Description']]
        if (!is.na(desc)) {
          desc <- unlist(x = strsplit(x = desc, split = '\n'))
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
        }
        # Assemble the information
        if (all(is.na(desc))) {
          desc <- c(
            'Dataset' = dataset,
            'Version' = pkg[['Version']],
            'Summary' = pkg[['Title']]
          )
        } else {
          desc <- c(
            'Dataset' = dataset,
            'Version' = pkg[['Version']],
            'Summary' = pkg[['Title']],
            desc
          )
        }
        return(desc)
      }
    )
    # Pad missing metadata with NAs
    manifest.names <- unique(x = unlist(
      x = lapply(X = avail.pkgs, FUN = names),
      use.names = FALSE
    ))
    for (pkg in names(x = avail.pkgs)) {
      for (col in manifest.names) {
        avail.pkgs[[pkg]][[col]] <- avail.pkgs[[pkg]][[col]] %||% NA
      }
    }
    # Convert each entry to a dataframe and bind everything together
    # avail.pkgs <- as.data.frame(t(avail.pkgs))
    avail.pkgs <- do.call(
      what = "rbind",
      args = sapply(X = avail.pkgs, FUN = as.data.frame, simplify = FALSE)
    )
    # Coerce version information to package_version
    avail.pkgs$Version <- package_version(x = avail.pkgs$Version)
  } else if (pkg.env$source == 'appdir') {
    # Read cached manifest
    pkg.env$extdata.warn <- FALSE
    packageStartupMessage(
      "Using cached data manifest, last updated at ",
      file.info(cache.manifest)$ctime
    )
    avail.pkgs <- readRDS(file = cache.manifest)
  } else if (pkg.env$source == 'extdata') {
    # Read SeuratData-bundled manifest
    if (!pkg.env$extdata.warn) {
      warning(
        "Using SeuratData-bundled data manifest. ",
        "This may be out-of-date and not contain the latest datasets. ",
        "This warning will be shown once per session or until we can read from a remote or cached data manifest",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    pkg.env$extdata.warn <- TRUE
    avail.pkgs <- readRDS(file = system.file(
      'extdata/manifest.Rds',
      package = 'SeuratData',
      mustWork = TRUE
    ))
  }
  # Get dataset installation status
  avail.pkgs$Installed <- vapply(
    X = rownames(x = avail.pkgs),
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
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
  # Coerce version information to package_version
  # Allow NAs to become effectively NA_pacakge_version_
  avail.pkgs$InstalledVersion <- package_version(
    x = avail.pkgs$InstalledVersion,
    strict = FALSE
  )
  # # TODO: remove these when we allow loading of processed datasets
  # cols.remove <- c('default.dataset', 'other.datasets')
  # if (any(cols.remove %in% colnames(x = avail.pkgs))) {
  #   ds.index <- which(x = colnames(x = avail.pkgs) %in% cols.remove)
  #   avail.pkgs <- avail.pkgs[, -ds.index]
  # }
  pkg.env$manifest <- avail.pkgs
  # Cache the manifest
  if (getOption(x = 'SeuratData.manifest.cache', default = FALSE)) {
    if (!dir.exists(paths = dirname(path = cache.manifest))) {
      dir.create(path = dirname(path = cache.manifest), recursive = TRUE)
    }
    cached <- if (file.exists(cache.manifest)) {
      readRDS(file = cache.manifest)
    } else {
      NULL
    }
    if (!isTRUE(x = all.equal(target = pkg.env$manifest, current = cached))) {
      saveRDS(object = pkg.env$manifest, file = cache.manifest)
    }
  }
  return(invisible(x = NULL))
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
