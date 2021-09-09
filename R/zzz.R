#' @import zlibbioc
#' @import Rcpp

UT_check_is_package <- getFromNamespace("check_is_package", "usethis")
UT_check_uses_roxygen <- getFromNamespace("check_uses_roxygen", "usethis")

UT_use_directory <- getFromNamespace("use_directory", "usethis")
UT_use_git_ignore <- getFromNamespace("use_git_ignore", "usethis")

UT_use_dependency <- getFromNamespace("use_dependency", "usethis")
UT_use_src <- getFromNamespace("use_src", "usethis")

UT_roxygen_ns_append <- getFromNamespace("roxygen_ns_append", "usethis")

#' ompBAM: C++ header library for parallel sequential reading
#' of BAM files
#' 
#' Although various tools for multi-threaded BAM processing (samtools,
#'	sambamba) are available, currently there are none available to the 
#'	R/Bioconductor environment. Parallelism in R is achieved using BiocParallel 
#'	or related packages. Typical alternatives include processing multiple BAM 
#'	files, each using a single core. Although easy to set up, memory-intensive 
#'  applications limit such approaches. ompBAM is a header-only C++ library 
#'	based on OpenMP, allowing developers to implement OpenMP-based parallelism 
#'	for the sequential reading of BAM files. ompBAM makes it easy by handling 
#'	all the parallel file-access and bgzf decompression, allowing developers to 
#'	focus on multi-threaded handling of individual reads.
#'
#' @details
#' For a quick start guide and a full list of functionality, refer to
#' the vignette
#' \code{vignette("ompBAM-API-Docs", package = "ompBAM")}
#'
#' @author Alex Chit Hei Wong
#' 
#' @docType package
#' @name ompBAM-package
#' @aliases ompBAM-package
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

#' Returns the path of a test BAM file
#'
#' Two types of BAM files are available. "Unsorted" is an unsorted BAM file of
#' paired RNA-seq of THP1 cell culture from Green et al (GSE130011, GSM3729250),
#' truncated to the first 10k reads. "scRNAseq" is a single cell 10X chromium V2
#' BAM file (sorted). It is sample E3.5_Rep1 from Nowotschin et al 
#' (GSE123046 / GSM3494334), truncated to the first 50k reads.
#'
#' @param dataset Returns the path to either the "Unsorted" or "scRNAseq" BAM.
#' @return A file path to the specified BAM file.
#' @examples
#' example_BAM("Unsorted")
#' @export
example_BAM <- function(dataset = c("Unsorted", "scRNAseq")) {
    dataset = match.arg(dataset)
    if(dataset == "Unsorted") return(
        system.file(file.path('extdata', 'THP1_ND_1.bam'), package = 'ompBAM'))
    if(dataset == "scRNAseq") return(
        system.file(file.path('extdata', 'MGE_E35_1.bam'), package = 'ompBAM'))
    stop("dataset must be one of 'Unsorted' or 'scRNAseq'")
}

#' Sets up the package in the given directory to use ompBAM
#'
#' This function creates a new package in the given directory (unless one 
#' already exists), then creates the necessary changes to run
#' ompBAM. It sets up the requisite dependencies and 'Make' files
#' so that the source code will successfully compile and link to the appropriate
#' external libraries (OpenMP, zlib) required for ompBAM.
#' @param path The path to the desired directory in which to set up a new
#'   R package. The directory name is assumed to be the package name.
#' @return None
#' @examples
#' path <- file.path(tempdir(), "myPkgName")
#' use_ompBAM(path)
#'
#' @export
use_ompBAM <- function(path = ".") {
    require("usethis")
    proj_path = .check_ompBAM_path(path)
    proj_name = basename(proj_path)   
    
    usethis::create_package(proj_path, open = FALSE)
    usethis::proj_set(proj_path, force = TRUE)
    
    makevars = file.path(proj_path, "src", "Makevars")
    makevars_win = file.path(proj_path, "src", "Makevars.win")
    example_code = file.path(proj_path, "src/ompBAM_example.cpp")
    R_import_file = file.path(proj_path, "R/ompBAM_imports.R")
    UT_use_directory("src")
    UT_use_git_ignore(c("*.o", "*.so", "*.dll"), "src")

    write(paste0("#' @useDynLib ", proj_name, ", .registration = TRUE"),
        R_import_file, sep="\n")
    write("#' @import ompBAM", R_import_file, append = TRUE, sep="\n")
    write("NULL", R_import_file, append = TRUE, sep="\n")
    ui_done(paste("Created", R_import_file, "with roxygen tags"))

    file.copy(system.file(
        file.path('examples', 'ompBAMExample', "src", "Makevars"), 
            package = 'ompBAM'),
        makevars)
    file.copy(system.file(
        file.path('examples', 'ompBAMExample', "src", "Makevars.win"), 
            package = 'ompBAM'),
        makevars_win)
    ui_done("Created src/Makevars and src/Makevars.win")

    file.copy(system.file(
        file.path('extdata', 'ompBAM_example.cpp'), package = 'ompBAM'),
        example_code)
    ui_done("Created src/ompBAM_example.cpp with idxstats_pbam() function")
    
    UT_use_dependency("ompBAM", "Imports")
    UT_use_dependency("Rcpp", "LinkingTo")
    UT_use_dependency("zlibbioc", "LinkingTo")
    UT_use_dependency("ompBAM", "LinkingTo")
    
    end_msg = paste(proj_name, "successfully created in", proj_path)
    ui_done(end_msg)
    
    end_msg2 = paste("Please run roxygen2 using devtools::document() before",
        "building the package.")
    message(end_msg2)
    
    
}

.check_ompBAM_path = function(path) {
    if(!dir.exists(dirname(path))) {
        errormsg = paste(dirname(path), "needs to exist")
        stop(errormsg, call. = FALSE)
    }
    if(dir.exists(path)) {
        proj_path = normalizePath(path)
    } else {
        proj_path = file.path(normalizePath(dirname(path)), basename(path))
    }
    makevars = file.path(proj_path, "src", "Makevars")
    makevars_win = file.path(proj_path, "src", "Makevars.win")
    if(file.exists(makevars) | file.exists(makevars_win)) {
        fail_make_msg = paste(
            "src/Makevars or src/Makevars.win already exist.",
            "use_ompBAM not run()"
        )
        stop(fail_make_msg, call. = FALSE)
    }
    example_code = file.path(proj_path, "src/ompBAM_example.cpp")
    if(file.exists(example_code)) {
        fail_src_msg = paste(
            example_code, "already exist.",
            "use_ompBAM not run()"
        )
        stop(fail_src_msg, call. = FALSE)
    }
    R_import_file = file.path(proj_path, "R/ompBAM_imports.R")
    if(file.exists(R_import_file)) {
        fail_src_msg = paste(
            R_import_file, "already exist.",
            "Has use_ompBAM() already been run?"
        )
        stop(fail_src_msg, call. = FALSE)
    }
    return(proj_path)
}