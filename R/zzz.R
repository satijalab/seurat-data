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

#' Index a saved \code{Seurat} object
#'
#' @inheritParams LoadObject
#'
#' @return ...
#'
#' @importFrom tools file_path_sans_ext
#'
#' @seealso \code{\link{SaveObject}} \code{\link{LoadObject}}
#'
#' @keywords internal
#'
IndexObject <- function(directory) {
  if (!dir.exists(paths = directory)) {
    stop("Cannot find directory ", directory, call. = FALSE)
  }
  # Get a file listing of the saved object directory
  files <- list.files(path = directory, full.names = TRUE, recursive = FALSE)
  files <- grep(pattern = '\\.Rds$', x = files, value = TRUE)
  if (length(x = files) < 1) {
    stop("No Rds files found in ", directory, call. = FALSE)
  }
  # Check Assays
  assay.dir <- file.path(directory, 'assays')
  if (!dir.exists(paths = assay.dir)) {
    stop("Cannot find assay directory", call. = FALSE)
  }
  assays <- list.dirs(path = assay.dir, full.names = FALSE, recursive = FALSE)
  assays <- sapply(
    X = assays,
    FUN = function(x) {
      return(grep(
        pattern = '\\.Rds$',
        x = list.files(
          path = file.path(assay.dir, x),
          full.names = TRUE,
          recursive = FALSE
        ),
        value = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  assays <- sapply(
    X = assays,
    FUN = function(x) {
      names(x = x) <- basename(path = file_path_sans_ext(x = x))
      return(x)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  # Check DimReducs
  reduc.dir <- file.path(directory, 'reductions')
  reducs <- list.dirs(path = reduc.dir, full.names = FALSE, recursive = FALSE)
  reducs <- sapply(
    X = reducs,
    FUN = function(x) {
      return(grep(
        pattern = '\\.Rds$',
        x = list.files(
          path = file.path(reduc.dir, x),
          full.names = TRUE,
          recursive = FALSE
        ),
        value = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  # Check Graphs
  # Check Misc
  # Check Tools
  return(invisible(x = NULL))
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
#' @keywords internal
#'
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

#' Load a saved \code{Seurat} object
#'
#' @param directory Path to the directory the dataset is stored
#' @param type How to load the \code{Seurat} object; choose from
#' \describe{
#'   \item{info}{Information about the object and what's stored in it}
#'   \item{raw}{The raw form of the dataset, no other options are evaluated}
#'   \item{processed}{...}
#' }
#' @param assays ...
#' @param dimreducs ...
#' @param graphs ...
#'
#' @return ...
#'
#' @seealso \code{\link{SaveObject}} \code{\link{IndexObject}}
#'
#' @keywords internal
#'
LoadObject <- function(
  directory,
  type = c('info', 'processed'),
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL
) {
  type <- match.arg(arg = type, choices = c('info', 'processed'))
  index <- IndexObject(directory = directory)
  if (type == 'info') {
    return(index)
  }
  return(invisible(x = NULL))
}

#' Find dataset packages based on name
#'
#' @param ds Names of datasets
#'
#' @return A vector of package names based on dataset names
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

#' Save a \code{Seurat} object into component Rds files
#'
#' @param object A \code{Seurat} object
#' @param directory Directory to save the object to
#' @param verbose Show progress messages
#'
#' @return Invisibly return the directory that the \code{Seurat} object was stored in
#'
#' @importFrom Seurat Project
#' @importFrom methods slotNames slot
#'
#' @seealso \code{\link{IndexObject}} \code{\link{LoadObject}}
#'
#' @keywords internal
#'
SaveObject <- function(
  object,
  directory = file.path(getwd(), Project(object = object)),
  verbose = TRUE
) {
  # Ensure the output directory exists
  if (!dir.exists(paths = directory)) {
    dir.create(path = directory, recursive = TRUE)
  }
  if (verbose) {
    message("Saving Seurat object to ", directory)
  }
  # Get the slots of the object
  names <- slotNames(x = object)
  # Figure out which are lists (eg. assays, graphs, reductions)
  collections <- names(x = which(x = sapply(
    X = names,
    FUN = function(n) {
      x <- slot(object = object, name = n)
      return(is.list(x = x) && !is.data.frame(x = x))
    }
  )))
  # For every collection
  for (col in collections) {
    obj <- slot(object = object, name = col)
    # Figure out which slots are S4 objects (eg. Assay, DimReduc, dgCMatrix)
    s4 <- vapply(
      X = obj,
      FUN = function(x) {
        return(isS4(x) && !inherits(x = x, what = c('Graph', 'JackStrawData', 'SeuratCommand')))
      },
      FUN.VALUE = logical(length = 1L)
    )
    if (length(x = s4) > 0) {
      # Determine if all values of obj are S4 objects (eg. list of Assay objects)
      if (all(s4)) {
        for (i in 1:length(x = obj)) {
          # Save each complex object in its own directory
          obj.name <- names(x = obj)[i]
          if (verbose) {
            message("Saving ", obj.name, " from ", col)
          }
          SaveObject(
            object = obj[[i]],
            directory = file.path(directory, col, obj.name),
            verbose = verbose
          )
        }
      } else {
        # Not a collection of S4 objects, just save the list (eg. misc, tools)
        if (length(x = obj) > 0) {
          if (verbose) {
            message("Saving ", col)
          }
          saveRDS(
            object = obj,
            file = file.path(directory, paste0(col, '.Rds'))
          )
        }
      }
    } else {
      # Save collections without S4 objects
      if (length(x = obj) > 0) {
        if (verbose) {
          message("Saving ", col)
        }
        saveRDS(
          object = obj,
          file = file.path(directory, paste0(col, '.Rds'))
        )
      }
    }
  }
  # Save everything that's not a list
  names <- setdiff(x = names, y = collections)
  for (x in names) {
    obj <- slot(object = object, name = x)
    if (!is.null(x = obj) && !IsMatrixEmpty(x = obj)) {
      if (verbose) {
        message("Saving ", x)
      }
      saveRDS(
        object = slot(object = object, name = x),
        file = file.path(directory, paste0(x, '.Rds'))
      )
    }
  }
  return(invisible(x = directory))
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
  # TODO: remove these when we allow loading of processed datasets
  ds.index <- which(x = colnames(x = avail.pkgs) %in% c('default.dataset', 'other.datasets'))
  avail.pkgs <- avail.pkgs[, -ds.index]
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
