#' ParaBAM: C++ header library for parallel sequential reading
#' of BAM files
#' 
#' ParaBAM provides a set of tools to facilitate developers to design
#' Rcpp / C++ based packages to read BAM files in a multi-threaded,
#' sequential approach. Currently, a multi-threaded solution is not
#' available to developers using Rhtslib.
#'
#' @details
#' For a quick start guide and a full list of functionality, refer to
#' the vignette
#' \code{vignette("ParaBAM-API-Docs", package = "ParaBAM")}
#'
#' @author Alex Chit Hei Wong
#' 
#' @docType package
#' @name NxtIRF-package
#' @aliases NxtIRF-package
#' @keywords package
#' @md
NULL

#' ParaBAM Example
#'
#' Installs the ParaBAMExample package included with the ParaBAM package
#'
#' @details
#' This function installs the ParaBAMExample package located inside
#' the 'examples' subfolder of ParaBAM.
#' @param pkg The name of the example package (default "ParaBAMExample")
#' @return None.
#' @examples
#' # The directory containing the source code is given by the path here
#' 
#' print(system.file(file.path('examples', "ParaBAMExample"), 
#'     package = 'ParaBAM'))
#' # Install the ParaBAMExample package
#' 
#' install_ParaBAM_example()
#' @export
install_ParaBAM_example <- function(pkg = "ParaBAMExample") {
    from <- system.file(file.path('examples', pkg), package = 'ParaBAM')
    if(!dir.exists(from)) {
        stop("Invalid example package name")
    }
    dir <- tempfile()
    dir.create(dir)
    file.copy(from, dir, recursive = TRUE)
    path <- file.path(dir, pkg)

    devtools::load_all(path, quiet = TRUE)
}