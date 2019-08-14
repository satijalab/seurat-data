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
#' \strong{Note}: Each reduction listed must be associated with an assay listed
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
#'
#' @return If \code{type} is info, information about the data contained within the
#' h5Seurat file. Otherwise, a \code{Seurat} object with the data requested
#'
#' @rdname LoadH5Seurat
#'
#' @seealso \code{\link{SaveH5Seurat}} \code{\link{IndexH5Seurat}}
#'
#' @keywords internal
#'
LoadH5Seurat <- function(file, ...) {
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
#' @rdname SaveH5Seurat
#'
#' @seealso \code{\link{LoadH5Seurat}} \code{\link{IndexH5Seurat}}
#'
#' @keywords internal
#'
SaveH5Seurat <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
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
#' @return A list with the data contained within \code{x}
#'
#' @aliases as.list
#'
#' @rdname as.list.H5Group
#' @method as.list H5Group
#'
#' @keywords internal
#'
as.list.H5Group <- function(x, ...) {
  return(invisible(x = NULL))
}

#' @rdname IndexH5Seurat
#' @method IndexH5Seurat character
#'
#' @keywords internal
#'
IndexH5Seurat.character <- function(file, ...) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to enable indexing h5Seurat files", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop("Cannot find h5Seurat file ", file, call = FALSE)
  }
  hfile <- hdf5r::H5File$new(filename = file, mode = 'r')
  on.exit(expr = hfile$close_all())
  return(IndexH5Seurat(file = hfile, ...))
}

#' @rdname IndexH5Seurat
#' @method IndexH5Seurat H5File
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
      slots <- c(
        'cell.embeddings',
        'feature.loadings',
        'feature.loadings.projected',
        'jackstraw'
      )
      check <- slots %in% names(x = file[['reductions']][[x]])
      names(x = check) <- slots
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
  # Get commands
  # Get metadata
  # Get misc
  # Get tools
  # Return the index
  return(structure(.Data = index, class = c('h5SI', 'list')))
}

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat character
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
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to enable loading Seurat objects from h5Seurat files", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop("Cannot find h5Seurat file ", file, call. = FALSE)
  }
  hfile <- hdf5r::H5File$new(filename = file, mode = 'r')
  on.exit(expr = hfile$close_all())
  return(LoadH5Seurat(
    x = hfile,
    type = type,
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    verbose = verbose,
    ...
  ))
}

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat H5File
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
  type <- match.arg(arg = type, choices = c('info', 'object'))
  # index <- IndexObject(directory = directory)
  # if (type == 'info') {
  #   return(index)
  # }
  # Load Assays
  # Load DimReducs
  if (!(isFALSE(x = reductions) || is.na(x = reductions))) {
    ''
  } else if (verbose) {
    ''
  }
  # Load Graphs
  if (!(isFALSE(x = graphs) || is.na(x = graphs))) {
    ''
  } else if (verbose) {
    ''
  }
  # Load SeuratCommands
  # Load meta.data
  # Load misc
  # Load tools
  return(invisible(x = NULL))
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
#'
#' @keywords internal
#'
print.h5SI <- function(x, ...) {
  catn <- function(...) {
    cat(..., '\n', sep = '')
  }
  symbols <- c(red(symbol$cross), green(symbol$tick))
  # Show assay information
  catn("Assays:")
  assays <- paste0(' ', names(x = x$assays), ': ')
  assays <- col_align(text = assays, width = max(nchar(x = assays)))
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
  return(invisible(x = x))
}

#' @importFrom Seurat as.Seurat Project
#'
#' @rdname SaveH5Seurat
#' @method SaveH5Seurat default
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
#' @importFrom Seurat Project DefaultAssay GetAssayData Key VariableFeatures Misc
#' Embeddings Loadings Stdev JS Command Tool
#'
#' @rdname SaveH5Seurat
#' @method SaveH5Seurat Seurat
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
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to enable saving Seurat object to h5Seurat files", call. = FALSE)
  }
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
  # Add Assays
  assays <- Filter(
    f = function(x) {
      return(inherits(x = object[[x]], what = 'Assay'))
    },
    x = names(x = object)
  )
  assay.group <- hfile$create_group(name = 'assays')
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
          dat.grp[['indices']] <- slot(object = dat, name = 'i') - 1
          dat.grp[['indptr']] <- slot(object = dat, name = 'p')
          dat.grp[['data']] <- slot(object = dat, name = 'x')
        } else {
          warning("wrong data")
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
          message("Adding miscellaneous data ", i, "for ", assay)
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
  reduc.group <- hfile$create_group(name = 'reductions')
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
          obj.group[[i]] <- Loadings(object = obj, projected = projected)
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
          message("Adding miscellaneous data ", i$name)
        }
        misc.group[[i$name]] <- i$value
      }
    }
  } else if (verbose) {
    message("No dimensional reduction information present")
  }
  # Add Graphs
  graph.group <- hfile$create_group(name = 'graphs')
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
  cmd.group <- hfile$create_group(name = 'commands')
  for (i in Command(object = object)) {
    i.group <- cmd.group$create_group(name = i)
    if (verbose) {
      message("Adding command ", i)
    }
    WriteH5List(
      x = Filter(
        f = Negate(f = is.function),
        x = as.list(x = Command(object = object, command = i))
      ),
      name = i,
      hgroup = i.group
    )
  }
  # Add meta.data
  hfile[['meta.data']] <- object[[]]
  # Add misc
  misc.group <- hfile$create_group(name = 'misc')
  for (m in Enumerate(x = Misc(object = object))) {
    if (verbose) {
      message("Adding miscellaneous data ", m$name)
    }
    WriteH5List(x = m$value, name = m$name, hgroup = misc.group)
  }
  # Add tools
  tool.group <- hfile$create_group(name = 'tools')
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

#' Write a list as a series of HDF5 groups and datasets
#'
#' @param x A list
#' @param name Name to save dataset as
#' @param hgroup An HDF5 file or group (\code{H5File} or \code{H5Group} objects from hdf5r)
#'
#' @return Invisibly returns \code{NULL}
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
      } else {
        xgroup[[i$name]] <- i$value
      }
    }
  } else {
    hgroup[[name]] <- x
  }
  return(invisible(x = NULL))
}
