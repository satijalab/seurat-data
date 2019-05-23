manifest <- vector(mode = 'character')

attached <- vector(mode = 'character')

.onLoad <- function(libname, pkgname) {
  UpdateManifest()
  for (pkg in manifest) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      message("Attaching ", pkg)
      try(expr = attachNamespace(ns = pkg), silent = TRUE)
      attached <<- c(attached, pkg)
    }
  }
}

.onUnload <- function(libpath) {
  for (pkg in attached) {
    message("Detaching ", pkg)
    unloadNamespace(ns = pkg)
  }
}