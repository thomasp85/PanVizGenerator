.onAttach <- function(libname, pkgname) {
    if (!hasGO()) {
        packageStartupMessage("No gene ontology file detected. Run 'getGO()' to download latest version")
    }
}