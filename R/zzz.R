#' ompBAM: C++ header library for parallel sequential reading
#' of BAM files
#' 
#' ompBAM provides a set of tools to facilitate developers to design
#' Rcpp / C++ based packages to read BAM files in a multi-threaded,
#' sequential approach. Currently, a multi-threaded solution is not
#' available to developers using Rhtslib.
#'
#' @details
#' For a quick start guide and a full list of functionality, refer to
#' the vignette
#' \code{vignette("ompBAM-API-Docs", package = "ompBAM")}
#'
#' @author Alex Chit Hei Wong
#' 
#' @docType package
#' @name NxtIRF-package
#' @aliases NxtIRF-package
#' @keywords package
#' @md
NULL

#' ompBAM Example
#'
#' Installs the ompBAMExample package included with the ompBAM package
#'
#' @details
#' This function installs the ompBAMExample package located inside
#' the 'examples' subfolder of ompBAM.
#' @param pkg The name of the example package (default "ompBAMExample")
#' @return None.
#' @examples
#' # The directory containing the source code is given by the path here
#' 
#' print(system.file(file.path('examples', "ompBAMExample"), 
#'     package = 'ompBAM'))
#' # Install the ompBAMExample package
#' 
#' install_ompBAM_example()
#' @export
install_ompBAM_example <- function(pkg = "ompBAMExample") {
    from <- system.file(file.path('examples', pkg), package = 'ompBAM')
    if(!dir.exists(from)) {
        stop("Invalid example package name")
    }
    dir <- tempfile()
    dir.create(dir)
    file.copy(from, dir, recursive = TRUE)
    path <- file.path(dir, pkg)

    devtools::load_all(path, quiet = TRUE)
}