repo.use <- 'http://satijalab04.nygenome.org'
pkg.env <- new.env()
pkg.env$manifest <- vector(mode = 'character')
pkg.env$attached <- vector(mode = 'character')

.onLoad <- function(libname, pkgname) {
  UpdateManifest()
  for (pkg in rownames(x = pkg.env$manifest)) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      message("Attaching ", pkg.env$manifest[pkg, 'Dataset', drop = TRUE])
      try(expr = attachNamespace(ns = pkg), silent = TRUE)
      pkg.env$attached <- c(pkg.env$attached, pkg)
    }
  }
}

.onUnload <- function(libpath) {
  for (pkg in pkg.env$attached) {
    message("Detaching ", pkg)
    unloadNamespace(ns = pkg)
  }
}
