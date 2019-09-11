#' @include zzz.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Index an H5Seurat object
#'
#' @param ... ...
#'
#' @return ...
#'
#' @rdname IndexH5Seurat
#'
#' @seealso \code{\link{LoadH5Seurat}} \code{\link{SaveH5Seurat}}
#'
#' @keywords internal
#'
IndexH5Seurat <- function(file, ...) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to enable indexing h5Seurat files", call. = FALSE)
  }
  UseMethod(generic = 'IndexH5Seurat', object = file)
}

#' Load a saved \code{Seurat} object from an h5Seurat file
#'
#' @param file Name of h5Seurat or connected h5Seurat file to load
#' @param assays One of:
#' \itemize{
#'   \item A character vector with names of assays
#'   \item A character vector with one or more of \code{counts}, \code{data},
#'   \code{scale.data} describing which slots of the \strong{default assay}
#'   \item A named list where each entry is either the name of an assay or a vector
#'   describing which slots (described above) to take from which assay
#'   \item \code{NULL} for all assays
#' }
#' @param reductions One of:
#' \itemize{
#'   \item A character vector with names of reductions
#'   \item \code{NULL} for all reductions
#'   \item \code{NA} or \code{FALSE} for no reductions
#' }
#' \strong{Note}: Only reductions associated with an assay loaded
#' in \code{assays}
#' @param graphs One of:
#' \itemize{
#'   \item A character vector with names of graphs
#'   \item \code{NULL} for all graphs
#'   \item \code{NA} or \code{FALSE} for no reductions
#' }
#' @param meta.data Load object metadata?
#' @param commands Load command information
#' @param misc Load miscellaneous data?
#' @param tools Load tool-specific information?
#' @param verbose Show progress updates
#'
#' @return If \code{type} is info, information about the data contained within the
#' h5Seurat file. Otherwise, a \code{Seurat} object with the data requested
#'
#' @export
#'
#' @seealso \code{\link{SaveH5Seurat}} \code{\link{IndexH5Seurat}}
#'
#' @keywords internal
#'
LoadH5Seurat <- function(file, ...) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to enable loading Seurat objects from h5Seurat files", call. = FALSE)
  }
  UseMethod(generic = 'LoadH5Seurat', object = file)
}

#' Save a \code{Seurat} object to a h5Seurat file
#'
#' @param object A \code{Seurat} object
#' @param filename Name for h5Seurat file
#' @param overwrite If an h5Seurat file of name \code{filename} exists, overwrite it?
#' @param verbose Show progress updates
#' @param ... Arguments passed to other methods
#'
#' @return Invisibly returns the name of the h5Seurat file
#'
#' @export
#'
#' @seealso \code{\link{LoadH5Seurat}} \code{\link{IndexH5Seurat}}
#'
#' @keywords internal
#'
SaveH5Seurat <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to enable saving Seurat object to h5Seurat files", call. = FALSE)
  }
  UseMethod(generic = 'SaveH5Seurat', object = object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read an HDF5 group into a ist
#'
#' @param x An \code{H5Group} object (from hdf5r)
#' @param ... Ignored
#'
#' @return A list with the data contained within \code{x}; if the HDF5 attribute
#' \code{s4class} is set and is a class, will return an object of class \code{s4class}
#' instead
#'
#' @aliases as.list
#'
#' @rdname as.list.H5Group
#' @method as.list H5Group
#' @export
#'
#' @keywords internal
#'
as.list.H5Group <- function(x, ...) {
  objs <- names(x = x)
  to.return <- vector(mode = 'list', length = length(x = objs))
  names(x = to.return) <- objs
  for (i in objs) {
    to.return[[i]] <- if (inherits(x = x[[i]], what = 'H5Group')) {
      as.list(x = x[[i]])
    } else if (length(x = x[[i]]$dims) > 1) {
      x[[i]][, ]
    } else {
      x[[i]][]
    }
  }
  if (!is.null(x = hdf5r::h5attributes(x = x)$s4class)) {
    to.return <- tryCatch(
      expr = do.call(
        what = 'new',
        args = c(Class = hdf5r::h5attr(x = x, which = 's4class'), to.return)
      ),
      error = function(...) {
        return(to.return)
      }
    )
  }
  return(to.return)
}

#' @rdname IndexH5Seurat
#' @method IndexH5Seurat character
#' @export
#'
#' @keywords internal
#'
IndexH5Seurat.character <- function(file, ...) {
  if (!file.exists(file)) {
    stop("Cannot find h5Seurat file ", file, call = FALSE)
  }
  hfile <- hdf5r::H5File$new(filename = file, mode = 'r')
  on.exit(expr = hfile$close_all())
  return(IndexH5Seurat(file = hfile, ...))
}

#' @rdname IndexH5Seurat
#' @method IndexH5Seurat H5File
#' @export
#'
#' @keywords internal
#'
IndexH5Seurat.H5File <- function(file, ...) {
  index <- list()
  # Get assay information
  assays <- names(x = file[['assays']])
  assays <- sapply(
    X = assays,
    FUN = function(x) {
      slots <- c('counts', 'data', 'scale.data')
      check <- slots %in% names(x = file[['assays']][[x]])
      names(x = check) <- slots
      check[['scale.data']] <- check[['scale.data']] && file[['assays']][[x]]$exists(name = 'scaled.features')
      return(check)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (!all(vapply(X = assays, FUN = '[[', FUN.VALUE = logical(length = 1L), 'data'))) {
    stop("All assays need a 'data' slot", call. = FALSE)
  }
  index[['assays']] <- assays
  # Get dimensional reduction information
  reductions <- names(x = file[['reductions']])
  reductions <- sapply(
    X = reductions,
    FUN = function(x) {
      reduc <- file[['reductions']][[x]]
      slots <- c(
        'cell.embeddings',
        'feature.loadings',
        'feature.loadings.projected',
        'jackstraw'
      )
      check <- slots %in% names(x = reduc)
      names(x = check) <- slots
      check[['feature.loadings']] <- check[['feature.loadings']] && reduc$exists(name = 'features')
      check[['feature.loadings.projected']] <- check[['feature.loadings.projected']] && reduc$exists(name = 'projected.features')
      check <- c(
        'assay' = hdf5r::h5attr(
          x = file[['reductions']][[x]],
          which = 'active.assay'
        ),
        as.list(x = check)
      )
      return(check)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  index[['reductions']] <- reductions
  # Get graph information
  graphs <- names(x = file[['graphs']])
  graph.assays <- sapply(X = strsplit(x = graphs, split = '_'), FUN = '[[', 1)
  graphs <- sapply(
    X = unique(x = graph.assays),
    FUN = function(x) {
      return(grep(pattern = paste0('^', x, '_'), x = graphs, value = TRUE))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  index[['graphs']] <- graphs
  # Get active.assay
  index[['active.assay']] <- hdf5r::h5attr(x = file, which = 'active.assay')
  # Get commands
  # Get metadata
  # Get misc
  # Get tools
  # Return the index
  return(ValidateH5SI(x = structure(.Data = index, class = 'h5SI')))
}

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat character
#' @export
#'
#' @keywords internal
#'
LoadH5Seurat.character <- function(
  file,
  type = c('info', 'object'),
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  verbose = TRUE,
  ...
) {
  if (!file.exists(file)) {
    stop("Cannot find h5Seurat file ", file, call. = FALSE)
  }
  hfile <- hdf5r::H5File$new(filename = file, mode = 'r')
  on.exit(expr = hfile$close_all())
  return(LoadH5Seurat(
    file = hfile,
    type = type,
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    verbose = verbose,
    ...
  ))
}

#' @importFrom methods new slot<-
#' @importFrom Seurat as.sparse CreateAssayObject Key<- GetAssayData SetAssayData
#' AddMetaData VariableFeatures<- CreateDimReducObject as.Graph Idents<-
#' Project<- Project Misc<-
#'
#' @rdname LoadH5Seurat
#' @method LoadH5Seurat H5File
#' @export
#'
#' @keywords internal
#'
LoadH5Seurat.H5File <- function(
  file,
  type = c('info', 'object'),
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  verbose = TRUE,
  ...
) {
  type <- match.arg(arg = tolower(x = type), choices = c('info', 'object'))
  cells <- file[['cell.names']][]
  index <- IndexH5Seurat(file = file)
  if (type == 'info') {
    return(index)
  }
  # Load Assays
  assay.slots <- c('counts', 'data', 'scale.data')
  assay.msg <- 'Assay specification must include either the name of an assay or one or more assay slots'
  assays.all <- is.null(x = assays)
  if (is.null(x = assays)) {
    assays <- names(x = index$assays)
  }
  if (!is.null(x = names(x = assays))) {
    assays <- as.list(x = assays)
  }
  if (!is.list(x = assays)) {
    if (any(assays %in% names(x = index$assays)) && any(assays %in% assay.slots)) {
      stop("Ambiguous assays", call. = FALSE)
    } else if (any(assays %in% names(x = index$assays))) {
      assays <- assays[assays %in% names(x = index$assays)]
      assays <- sapply(
        X = assays,
        FUN = function(...) {
          return(assay.slots)
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
    } else {
      assays <- assays[assays %in% assay.slots]
      if (length(x = assays) < 1) {
        stop(assay.msg, call. = FALSE)
      }
      assays <- list(assays)
      assays <- rep_len(x = assays, length.out = length(x = index$assays))
      names(x = assays) <- names(x = index$assays)
    }
  } else {
    for (i in 1:length(x = assays)) {
      assay.name <- names(x = assays)[i] %||% names(x = index$assays)[i] %||% ''
      if (!assay.name %in% names(x = index$assays)) {
        if (assays[[i]][1] %in% names(x = index$assays)) {
          assay.name <- assays[[i]][1]
        } else if (any(assays[[i]] %in% assay.slots)) {
          assay.name <- hdf5r::h5attr(x = file, which = 'active.assay')
        }
      }
      if (nchar(x = assay.name) < 0 || !assay.name %in% names(x = index$assays)) {
        stop(assay.msg, call. = FALSE)
      }
      assay.content <- assays[[i]]
      if (assay.content[1] %in% names(x = index$assays)) {
        assay.content <- assay.slots
      } else {
        assay.content <- assay.content[assay.content %in% assay.slots]
        if (length(x = assay.content) < 1) {
          stop(assay.msg, call. = FALSE)
        }
      }
      assays[i] <- list(assay.content)
      names(x = assays)[i] <- assay.name
    }
  }
  assays.checked <- assays
  unique.assays <- unique(x = names(x = assays.checked))
  assays <- vector(mode = 'list', length = length(x = unique.assays))
  names(x = assays) <- unique.assays
  for (i in unique.assays) {
    assays.use <- which(x = names(x = assays.checked) == i)
    slots.use <- unique(x = unlist(x = assays.checked[assays.use], use.names = FALSE))
    slots.use <- slots.use[match(x = names(x = index$assays[[i]]), table = slots.use)]
    slots.use <- na.omit(object = slots.use[index$assays[[i]]])
    assays[[i]] <- slots.use
  }
  for (i in Enumerate(x = assays)) {
    if (!any(c('counts', 'data') %in% i$value)) {
      stop("All assays must have either a counts or data slot, missing for ", i$name, call. = FALSE)
    }
  }
  if (!index$active.assay %in% names(x = assays)) {
    warning(
      "Default assay not requested, using ",
      names(x = assays)[1],
      " as default instead",
      call. = FALSE,
      immediate. = TRUE
    )
    active.assay <- names(x = assays)[1]
  } else {
    active.assay <- index$active.assay
  }
  assays <- sapply(
    X = names(x = assays),
    FUN = function(x) {
      assay.group <- file[['assays']][[x]]
      features <- assay.group[['feature.names']][]
      # Add counts if not data, otherwise add data
      if ('counts' %in% assays[[x]] && !'data' %in% assays[[x]]) {
        if (verbose) {
          message("Initializing ", x, " with counts")
        }
        counts.data <- assay.group[['counts']]
        counts <- if (inherits(x = counts.data, what = 'H5Group')) {
          as.sparse(x = counts.data)
        } else {
          counts.data[, ]
        }
        rownames(x = counts.data) <- features
        colnames(x = counts.data) <- cells
        obj <- CreateAssayObject(
          counts = counts.data,
          min.cells = -1,
          min.features = -1
        )
      } else {
        if (verbose) {
          message("Initializing ", x, " with data")
        }
        data.data <- assay.group[['data']]
        data <- if (inherits(x = data.data, what = 'H5Group')) {
          as.sparse(x = data.data)
        } else {
          data.data[, ]
        }
        rownames(x = data) <- features
        colnames(x = data) <- cells
        obj <- CreateAssayObject(data = data)
      }
      Key(object = obj) <- hdf5r::h5attr(x = assay.group, which = 'key')
      # Add remaining slots
      for (slot in assays[[x]]) {
        if (IsMatrixEmpty(x = GetAssayData(object = obj, slot = slot))) {
          if (verbose) {
            message("Adding ", slot, " for ", x)
          }
          dat.data <- assay.group[[slot]]
          dat <- if (inherits(x = dat.data, what = 'H5Group')) {
            as.sparse(x = dat.data)
          } else {
            dat.data[, ]
          }
          rownames(x = dat) <- if (slot == 'scale.data') {
            assay.group[['scaled.features']][]
          } else {
            features
          }
          colnames(x = dat) <- cells
          obj <- SetAssayData(object = obj, slot = slot, new.data = dat)
        }
      }
      # Add meta features
      if (assay.group$exists(name = 'meta.features')) {
        if (verbose) {
          message("Adding feature-level metadata for ", x)
        }
        meta.data <- assay.group[['meta.features']][]
        rownames(x = meta.data) <- features
        obj <- AddMetaData(
          object = obj,
          metadata = meta.data
        )
      } else if (verbose) {
        message("No feature-level metadata for ", x)
      }
      # Add variable feature information
      if (assay.group$exists(name = 'variable.features')) {
        if (verbose) {
          message("Adding variable feature information for ", x)
        }
        VariableFeatures(object = obj) <- assay.group[['variable.features']][]
      } else if (verbose) {
        message("No variable features for ", x)
      }
      # Add miscellaneous information
      if (assay.group$exists(name = 'misc')) {
        slot(object = obj, name = 'misc') <- as.list(x = assay.group[['misc']])
      } else if (verbose) {
        message("No miscellaneous information for ", x)
      }
      return(obj)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  default.assays <- list(assays[[active.assay]])
  names(x = default.assays) <- active.assay
  object <- new(
    Class = 'Seurat',
    assays = default.assays,
    active.assay = active.assay,
    meta.data = data.frame(row.names = cells),
    version = package_version(x = hdf5r::h5attr(x = file, which = 'version'))
  )
  for (assay in names(x = assays)) {
    if (assay != active.assay) {
      object[[assay]] <- assays[[assay]]
    }
  }
  # Load DimReducs
  reducs.all <- is.null(x = reductions)
  if (is.null(x = reductions)) {
    reductions <- names(x = index$reductions)
  }
  if (!(isFALSE(x = reductions) || is.na(x = reductions))) {
    reductions <- Filter(
      f = function(reduc) {
        return(index$reductions[[reduc]]$assay %in% names(x = assays))
      },
      x = reductions
    )
    if (length(x = reductions) < 1) {
      warning(
        "None of the reductions specified are associated with a loaded assay",
        call. = FALSE,
        immediate. = TRUE
      )
    } else {
      reductions <- sapply(
        X = reductions,
        FUN = function(x) {
          reduc.group <- file[['reductions']][[x]]
          key <- hdf5r::h5attr(x = reduc.group, which = 'key')
          # Pull cell embeddings
          if (reduc.group$exists(name = 'cell.embeddings')) {
            if (verbose) {
              message("Pulling cell embeddings for ", x)
            }
            cell.embeddings <- reduc.group[['cell.embeddings']][, ]
            rownames(x = cell.embeddings) <- cells
            colnames(x = cell.embeddings) <- paste0(key, 1:ncol(x = cell.embeddings))
          } else {
            if (verbose) {
              message("No cell embeddings for ", x)
            }
            cell.embeddings <- new(Class = 'matrix')
          }
          # Pull feature loadings
          if (index$reductions[[x]]$feature.loadings) {
            if (verbose) {
              message("Pulling feature loadings for ", x)
            }
            feature.loadings <- reduc.group[['feature.loadings']][, ]
            rownames(x = feature.loadings) <- reduc.group[['features']][]
            colnames(x = feature.loadings) <- paste0(key, 1:ncol(x = feature.loadings))
          } else {
            if (verbose) {
              message("No feature loadings for ", x)
            }
            feature.loadings <- new(Class = 'matrix')
          }
          # Pull projected loadings
          if (index$reductions[[x]]$feature.loadings.projected) {
            if (verbose) {
              message("Pulling projected loadings for ", x)
            }
            projected.loadings <- reduc.group[['feature.loadings.projected']][, ]
            rownames(x = projected.loadings) <- reduc.group[['projected.features']][]
            colnames(x = projected.loadings) <- paste0(key, 1:ncol(x = projected.loadings))
          } else {
            if (verbose) {
              message("No projected loadings for ", x)
            }
            projected.loadings <- new(Class = 'matrix')
          }
          obj <- CreateDimReducObject(
            embeddings = cell.embeddings,
            loadings = feature.loadings,
            projected = projected.loadings,
            assay = index$reductions[[x]]$assay,
            stdev = if (reduc.group$exists(name = 'stdev')) {
              reduc.group[['stdev']][]
            } else {
              numeric()
            },
            key = key
          )
          # Add misc
          if (reduc.group$exists(name = 'misc')) {
            if (verbose) {
              message("Pulling miscellaneous information for ", x)
            }
            slot(object = obj, name = 'misc') <- as.list(x = reduc.group[['misc']])
          } else if (verbose) {
            message("No miscellaneous information for ", x)
          }
          # Add JackStraw
          if (index$reductions[[x]]$jackstraw) {
            if (verbose) {
              message("Pulling JackStraw information for ", x)
            }
            js <- new(Class = 'JackStrawData')
            for (slot in names(x = reduc.group[['jackstraw']])) {
              JS(object = js, slot = slot) <- as.matrix(x = reduc.group[['jackstraw']][[slot]][, ])
            }
            JS(object = obj) <- js
          } else if (verbose) {
            message("No JackStraw information for ", x)
          }
          return(obj)
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
      for (reduc in names(x = reductions)) {
        object[[reduc]] <- reductions[[reduc]]
      }
    }
  } else if (verbose) {
    message("No dimensional reduction information requested")
  }
  # Load Graphs
  graphs.all <- is.null(x = graphs)
  if (is.null(x = graphs)) {
    graphs.assay <- vapply(
      X = names(x = index$graphs),
      FUN = function(x) {
        return(any(grepl(pattern = x, x = names(x = assays), ignore.case = TRUE)))
      },
      FUN.VALUE = logical(length = 1L)
    )
    graphs <- unlist(x = index$graphs[graphs.assay])
    if (length(x = graphs) == 0) {
      graphs <- FALSE
    }
  }
  if (!(isFALSE(x = graphs) || is.na(x = graphs))) {
    for (graph in graphs) {
      if (verbose) {
        message("Adding graph ", graph)
      }
      graph.matrix <- as.sparse(x = file[['graphs']][[graph]])
      rownames(x = graph.matrix) <- colnames(x = graph.matrix) <- cells
      object[[graph]] <- as.Graph(x = graph.matrix)
    }
  } else if (verbose) {
    message("No nearest neighbor graphs requested")
  }
  # Load SeuratCommands
  if (all(assays.all, reducs.all, graphs.all)) {
    cmds <- vector(mode = 'list', length = length(x = names(x = file[['commands']])))
    names(x = cmds) <- names(x = file[['commands']])
    for (cmd in names(x = file[['commands']])) {
      if (verbose) {
        message("Loading command information for ", cmd)
      }
      cmds[[cmd]] <- new(
        Class = 'SeuratCommand',
        name = hdf5r::h5attr(x = file[['commands']][[cmd]], which = 'name'),
        time.stamp = as.POSIXct(x = hdf5r::h5attr(
          x = file[['commands']][[cmd]],
          which = 'time.stamp'
        )),
        call.string = hdf5r::h5attr(
          x = file[['commands']][[cmd]],
          which = 'call.string'
        ),
        params = as.list(x = file[['commands']][[cmd]])
      )
    }
    slot(object = object, name = 'commands') <- cmds
  }
  # Load meta.data
  meta.data <- file[['meta.data']][]
  rownames(x = meta.data) <- cells
  object <- AddMetaData(object = object, metadata = meta.data)
  # Set identity class and Project
  if (!is.null(x = hdf5r::h5attributes(x = file)$project)) {
    Project(object = object) <- hdf5r::h5attr(x = file, which = 'project')
  }
  if (file$exists(name = 'active.ident')) {
    idents <- file[['active.ident']][]
    names(x = idents) <- cells
    Idents(object = object) <- idents
  } else {
    Idents(object = object) <- Project(object = object)
  }
  # Load misc
  for (x in names(x = file[['misc']])) {
    if (verbose) {
      message("Adding miscellaneous data ", x)
    }
    Misc(object = object, slot = x) <- as.list(x = file[['misc']][[x]])
  }
  # Load tools
  slot(object = object, name = 'tools') <- as.list(x = file[['tools']])
  return(object)
}

#' Print index information for an h5Seurat Index (h5SI)
#'
#' @inheritParams base::print
#'
#' @return Invisibly returns \code{x}
#'
#' @importFrom cli symbol
#' @importFrom crayon red green col_align
#'
#' @aliases print
#'
#' @rdname print.h5SI
#' @method print h5SI
#' @export
#'
#' @keywords internal
#'
print.h5SI <- function(x, ...) {
  catn <- function(...) {
    cat(..., '\n', sep = '')
  }
  symbols <- c(red(symbol$cross), green(symbol$tick))
  x <- ValidateH5SI(x = x)
  # Show assay information
  catn("Assays:")
  assays <- names(x = x$assays)
  assays <- paste0(
    ' ',
    ifelse(test = assays == x$active.assay, yes = '*', no = ''),
    assays,
    ': '
  )
  assays <- col_align(
    text = assays,
    width = max(nchar(x = assays)),
    align = 'right'
  )
  catn(
    MakeSpace(n = max(nchar(x = assays))),
    col_align(
      text = c('counts', 'data', 'scale.data'),
      width = nchar(x = 'scale.data') + 1,
      align = 'center'
    )
  )
  for (i in 1:length(x = x$assays)) {
    catn(
      assays[i],
      col_align(
        text = symbols[x$assays[[i]] + 1],
        width = nchar(x = 'scale.data') + 1,
        align = 'center'
      )
    )
  }
  # Show reduction informaiton
  reduc.header <- c(
    'Embeddings' = 'cell.embeddings',
    'Loadings' = 'feature.loadings',
    'Projected' = 'feature.loadings.projected',
    'JackStraw' = 'jackstraw'
  )
  grouped.reductions <- split(
    x = x$reductions,
    f = unique(x = sapply(X = x$reductions, FUN = '[[', 'assay'))
  )
  for (assay in names(x = grouped.reductions)) {
    catn("Reductions for ", assay, ":")
    reduc.use <- grouped.reductions[[assay]]
    reductions <- paste0(' ', names(x = reduc.use), ': ')
    reductions <- col_align(text = reductions, width = max(nchar(x = reductions)))
    catn(
      MakeSpace(n = max(nchar(x = reductions))),
      col_align(
        text = names(x = reduc.header),
        width = max(nchar(x = names(x = reduc.header))) + 1,
        align = 'center'
      )
    )
    for (i in 1:length(x = grouped.reductions)) {
      catn(
        reductions[i],
        col_align(
          text = symbols[unlist(x = reduc.use[[i]][reduc.header]) + 1],
          width = max(nchar(x = names(x = reduc.header))) + 1,
          align = 'center'
        )
      )
    }
  }
  # Show graph information
  for (assay in names(x = x$graphs)) {
    catn("Graphs for ", assay, ":")
    catn(' ', paste(x$graphs[[assay]], collapse = ', '))
  }
  return(invisible(x = x))
}

#' @importFrom Seurat as.Seurat Project
#'
#' @rdname SaveH5Seurat
#' @method SaveH5Seurat default
#' @export
#'
#' @keywords internal
#'
SaveH5Seurat.default <- function(
  object,
  filename,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  tryCatch(
    expr = object <- as.Seurat(x = object, verbose = verbose, ...),
    error = function(...) {
      stop(
        "Unable to coerce an object of class ",
        class(x = object),
        " to a Seurat object",
        call. = FALSE
      )
    }
  )
  if (missing(x = filename)) {
    filename <- paste0(Project(object = object), '.h5Seurat')
  }
  return(invisible(x = SaveH5Seurat(
    object = object,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose
  )))
}

#' @importFrom tools file_ext
#' @importFrom methods slot slotNames
#' @importFrom Seurat Project DefaultAssay Idents GetAssayData Key
#' VariableFeatures Misc Embeddings Loadings Stdev JS Command Tool
#'
#' @rdname SaveH5Seurat
#' @method SaveH5Seurat Seurat
#' @export
#'
#' @keywords internal
#'
SaveH5Seurat.Seurat <- function(
  object,
  filename = paste0(Project(object = object), '.h5Seurat'),
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  if (file_ext(x = filename) != 'h5Seurat') {
    filename <- paste0(filename, 'h5Seurat')
  }
  if (file.exists(filename)) {
    if (overwrite) {
      warning(
        "Overwriting previous file ",
        filename,
        call. = FALSE,
        immediate. = TRUE
      )
      file.remove(filename)
    } else {
      stop("H5Seurat file at ", filename, " exists", call. = FALSE)
    }
  }
  hfile <- hdf5r::H5File$new(filename = filename, mode = 'w-')
  on.exit(expr = hfile$close_all())
  # Add attributes for project, version, and default assay
  hdf5r::h5attr(x = hfile, which = 'project') <- Project(object = object)
  hdf5r::h5attr(x = hfile, which = 'version') <- as.character(x = slot(
    object = object,
    name = 'version'
  ))
  hdf5r::h5attr(x = hfile, which = 'active.assay') <- DefaultAssay(object = object)
  # Add meta.data, identity class, and cell names
  hfile[['meta.data']] <- object[[]]
  hfile[['active.ident']] <- Idents(object = object)
  hfile[['cell.names']] <- colnames(x = object)
  # Prepare group information
  required.groups <- c('assays', 'reductions', 'graphs', 'commands', 'tools', 'misc')
  for (group in required.groups) {
    hfile$create_group(name = group)
  }
  # Add Assays
  assays <- Filter(
    f = function(x) {
      return(inherits(x = object[[x]], what = 'Assay'))
    },
    x = names(x = object)
  )
  assay.group <- hfile[['assays']]
  for (assay in assays) {
    obj <- object[[assay]]
    obj.group <- assay.group$create_group(name = assay)
    # Add counts/data/scale.data
    for (i in c('counts', 'data', 'scale.data')) {
      dat <- GetAssayData(object = obj, slot = i)
      if (!IsMatrixEmpty(x = dat)) {
        if (verbose) {
          message("Adding ", i, " for ", assay)
        }
        if (inherits(x = dat, what = 'matrix')) {
          obj.group[[i]] <- dat
        } else if (inherits(x = dat, what = 'dgCMatrix')) {
          dat.grp <- obj.group$create_group(name = i)
          dat.grp[['indices']] <- slot(object = dat, name = 'i')
          dat.grp[['indptr']] <- slot(object = dat, name = 'p')
          dat.grp[['data']] <- slot(object = dat, name = 'x')
        } else {
          warning("wrong data")
        }
        if (i == 'scale.data') {
          obj.group[['scaled.features']] <- rownames(x = dat)
        }
      } else if (verbose) {
        message("No data found in ", i, " for ", assay)
      }
    }
    # Add key
    hdf5r::h5attr(x = obj.group, which = 'key') <- Key(object = obj)
    # Add variable features
    if (length(x = VariableFeatures(object = obj)) > 0) {
      if (verbose) {
        message("Adding variable features for ", assay)
      }
      obj.group[['variable.features']] <- VariableFeatures(object = obj)
    } else if (verbose) {
      message("No variable features found for ", assay)
    }
    # Add meta.features
    if (ncol(x = obj[[]]) > 0) {
      if (verbose) {
        message("Adding meta.feature information for ", assay)
      }
      obj.group[['meta.features']] <- obj[[]]
    } else if (verbose) {
      message("No meta.feature information for ", assay)
    }
    # Add misc
    if (length(x = Misc(object = obj)) > 0) {
      misc.group <- obj.group$create_group(name = 'misc')
      for (i in names(x = Misc(object = obj))) {
        if (verbose) {
          message("Adding miscellaneous data: ", i, "for ", assay)
        }
        misc.group[[i]] <- Misc(object = obj, slot = i)
      }
    } else if (verbose) {
      message("No miscellaneous data found for ", assay)
    }
    # Add feature names
    obj.group[['feature.names']] <- rownames(x = obj)
  }
  # Add DimReducs
  reduc.group <- hfile[['reductions']]
  reductions <- Filter(
    f = function(x) {
      return(inherits(x = object[[x]], what = 'DimReduc'))
    },
    x = names(x = object)
  )
  if (length(x = reductions) > 0) {
    for (reduc in reductions) {
      obj <- object[[reduc]]
      obj.group <- reduc.group$create_group(name = reduc)
      # Add cell embeddings
      if (!IsMatrixEmpty(x = Embeddings(object = obj))) {
        if (verbose) {
          message("Adding cell embeddings for ", reduc)
        }
        obj.group[['cell.embeddings']] <- Embeddings(object = obj)
      } else if (verbose) {
        message("No cell embeddings for ", reduc)
      }
      # Add feature loadings
      for (i in c('feature.loadings', 'feature.loadings.projected')) {
        projected <- grepl(pattern = 'projected', x = i)
        type <- ifelse(test = projected, yes = 'projected loadings', no = 'loadings')
        if (!IsMatrixEmpty(x = Loadings(object = obj, projected = projected))) {
          if (verbose) {
            message("Adding ", type, " for ", reduc)
          }
          loadings <- Loadings(object = obj, projected = projected)
          obj.group[[i]] <- loadings
          features <- ifelse(
            test = projected,
            yes = 'projected.features',
            no = 'features'
          )
          obj.group[[features]] <- rownames(x = loadings)
        } else if (verbose) {
          message("No ", type, " for ", reduc)
        }
      }
      # Add default assay and key
      hdf5r::h5attr(x = obj.group, which = 'active.assay') <- DefaultAssay(object = obj)
      hdf5r::h5attr(x = obj.group, which = 'key') <- Key(object = obj)
      # Add standard deviations
      if (length(x = Stdev(object = obj)) > 0) {
        if (verbose) {
          message("Adding standard deviations for ", reduc)
        }
        obj.group[['stdev']] <- Stdev(object = obj)
      } else if (verbose) {
        message("No standard deviations for ", reduc)
      }
      # Add JackStraw
      if (as.logical(x = JS(object = obj))) {
        if (verbose) {
          message("Adding JackStraw information for ", reduc)
        }
        js.group <- obj.group$create_group(name = 'jackstraw')
        for (x in slotNames(x = JS(object = obj))) {
          js.group[[x]] <- JS(object = obj, slot = x)
        }
      } else if (verbose) {
        message("No JackStraw data for ", reduc)
      }
      # Add misc
      # TODO: Add Misc method for DimReducs
      misc.group <- obj.group$create_group(name = 'misc')
      for (i in Enumerate(x = slot(object = obj, name = 'misc'))) {
        if (verbose) {
          message("Adding miscellaneous data: ", i$name)
        }
        misc.group[[i$name]] <- i$value
      }
    }
  } else if (verbose) {
    message("No dimensional reduction information present")
  }
  # Add Graphs
  graph.group <- hfile[['graphs']]
  graphs <- Filter(
    f = function(x) {
      return(inherits(x = object[[x]], what = 'Graph'))
    },
    x = names(x = object)
  )
  if (length(x = graphs) > 0) {
    for (graph in graphs) {
      if (verbose) {
        message("Saving graph ", graph)
      }
      obj.group <- graph.group$create_group(name = graph)
      obj.group[['indices']] <- slot(object = object[[graph]], name = 'i')
      obj.group[['indptr']] <- slot(object = object[[graph]], name = 'p')
      obj.group[['data']] <- slot(object = object[[graph]], name = 'x')
    }
  } else if (verbose) {
    message("No nearest-neighbor graphs present")
  }
  # Add SeuratCommands
  cmd.group <- hfile[['commands']]
  for (i in Command(object = object)) {
    if (verbose) {
      message("Adding command ", i)
    }
    WriteH5List(
      x = Filter(
        f = Negate(f = is.function),
        x = as.list(x = Command(object = object, command = i))
      ),
      name = i,
      hgroup = cmd.group
    )
    cmd <- Command(object = object, command = i)
    for (slot in slotNames(x = cmd)) {
      if (slot == 'params') {
        next
      }
      hdf5r::h5attr(x = cmd.group[[i]], which = slot) <- as.character(x = slot(
        object = cmd,
        name = slot
      ))
    }
  }
  # Add misc
  misc.group <- hfile[['misc']]
  for (m in Enumerate(x = Misc(object = object))) {
    if (verbose) {
      message("Adding miscellaneous data: ", m$name)
    }
    WriteH5List(x = m$value, name = m$name, hgroup = misc.group)
  }
  # Add tools
  tool.group <- hfile[['tools']]
  for (x in Tool(object = object)) {
    if (verbose) {
      message("Writing data for tool ", x)
    }
    WriteH5List(x = Tool(object = object, slot = x), name = x, hgroup = tool.group)
  }
  return(invisible(x = filename))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Standard functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h5SI <- function(...) {
  object <- list(
    'assays' = list(),
    'reductions' = list(),
    'graphs' = list(),
    'active.assay' = vector(mode = 'character', length = 1L),
    'meta.data' = vector(mode = 'character'),
    ''
  )
  return(structure(.Data = object, class = 'h5SI'))
}

#' Validate an h5Seurat Index
#'
#' @param x An h5Seurat Index (h5SI)
#'
#' @return \code{x} if valid, otherwise a modified \code{x} to be valid
#'
#' @keywords internal
#'
ValidateH5SI <- function(x) {
  if (!inherits(x = x, what = 'h5SI')) {
    stop("'x' must be an h5Seurat index", call. = FALSE)
  }
  return(x)
}

#' Write a list as a series of HDF5 groups and datasets
#'
#' @param x A list
#' @param name Name to save dataset as
#' @param hgroup An HDF5 file or group (\code{H5File} or
#' \code{H5Group} objects from hdf5r)
#'
#' @return Invisibly returns \code{NULL}
#'
#' @importFrom methods slotNames slot
#'
#' @keywords internal
#'
WriteH5List <- function(x, name, hgroup) {
  if (!inherits(x = hgroup, what = c('H5File', 'H5Group'))) {
    stop("'hgroup' must be an H5File or H5Group object from hdf5r", call. = FALSE)
  }
  if (is.list(x = x) && !is.data.frame(x = x)) {
    xgroup <- hgroup$create_group(name = name)
    for (i in Enumerate(x = x)) {
      if (is.list(x = i$value) && !is.data.frame(x = i$value)) {
        WriteH5List(x = i$value, name = i$name, hgroup = xgroup)
      } else if (!is.null(x = i$value)) {
        xgroup[[i$name]] <- i$value
      }
    }
  } else if (isS4(x)) {
    xgroup <- hgroup$create_group(name = name)
    hdf5r::h5attr(x = xgroup, which = 's4class') <- class(x = x)[1]
    for (i in slotNames(x = x)) {
      obj <- slot(object = x, name = i)
      if (is.list(x = obj) && !is.data.frame(x = obj)) {
        WriteH5List(x = obj, name = i, hgroup = xgroup)
      } else if (!is.null(x = obj)) {
        xgroup[[i]] <- obj
      }
    }
  } else if (!is.null(x = x)) {
    hgroup[[name]] <- x
  }
  return(invisible(x = NULL))
}
