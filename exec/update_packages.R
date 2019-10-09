#!/usr/bin/env Rscript

# Check required packages
if (!requireNamespace('tools', quietly = TRUE)) {
  stop("Cannot find the tools package, exiting", call. = FALSE)
}

# Get the script name
ProgramName <- function() {
    prefix <- '--file='
    name <- sub(
      pattern = prefix,
      replacement = '',
      x = grep(
        pattern = paste0(prefix, "(.+)"),
        x = commandArgs(),
        value = TRUE
      )
    )
    return(name)
}

# Usage message
usage <- paste("Usage:", ProgramName(), "/path/to/repo")

# Parse and check arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(x = args) != 1 || tolower(x = args) %in% c('-h', '--help')) {
  cat(usage, '\n', file = stderr())
  quit(status = 1)
}
if (!dir.exists(paths = args)) {
  stop("Cannot find repo directory ", args, call. = FALSE)
}

# Get the contrib directories
dirs <- list.dirs(path = args)
dirs <- grep(pattern = 'contrib', x = dirs, value = TRUE)
if (length(x = dirs) < 1) {
  stop("Invalid repository structure in ", args, call. = FALSE)
}

# Write the PACKAGES manifests
for (d in dirs) {
  type <- if (grepl(pattern = 'macosx', x = d)) {
    'mac.binary'
  } else if (grepl(pattern = 'windows', x = d)) {
    'win.binary'
  } else {
    'source'
  }
  tools::write_PACKAGES(
    dir = d,
    fields = c('Description', 'Title', 'Date'),
    type = type,
    verbose = TRUE
  )
}
