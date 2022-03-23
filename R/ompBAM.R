#' @import zlibbioc
#' @import Rcpp
NULL

#' ompBAM: C++ header library for parallel sequential reading
#' of BAM files
#' 
#' Although various tools for multi-threaded BAM processing (samtools,
#'  sambamba) are available, currently there are none available to the 
#'  R/Bioconductor environment. Parallelism in R is achieved using BiocParallel 
#'  or related packages. Typical alternatives include processing multiple BAM 
#'  files, each using a single core. Although easy to set up, memory-intensive 
#'  applications limit such approaches. ompBAM is a header-only C++ library 
#'  based on OpenMP, allowing developers to implement OpenMP-based parallelism 
#'  for the sequential reading of BAM files. ompBAM makes it easy by handling 
#'  all the parallel file-access and bgzf decompression, allowing developers to 
#'  focus on multi-threaded handling of individual reads.
#'
#' @details
#' For a quick start guide and a full list of functionality, refer to
#' the vignette
#' \href{../doc/ompBAM-API-Docs.html}{\code{vignette("ompBAM-API-Docs", 
#' package = "ompBAM")}}
#'
#' @author Alex Chit Hei Wong
#' 
#' @docType package
#' @name ompBAM-package
#' @aliases ompBAM-package
#' @keywords package
#' @md
NULL

#' Sets up the package in the given directory to use ompBAM
#'
#' This function creates a new package in the given directory (unless one 
#' already exists).
#' It then sets up the requisite dependencies and 'Make' and 'configure' files
#' so that the source code will successfully compile and link to the appropriate
#' external libraries (OpenMP, zlib) required for ompBAM. All 3 major platforms
#' (Linux, MacOS and Windows) are supported, although MacOS users must first
#' install the required OpenMP libraries first. See details below.
#' @details OpenMP is natively supported on Linux and Windows but 
#'   not on MacOS. To compile ompBAM_based
#'   packages with OpenMP on MacOS, users should install the requisite OpenMP
#'   libraries. An easy way to do this would be to run: `brew install libomp`.
#'    For details, refer to 
#'   [this guide](https://mac.r-project.org/openmp/)
#' @param path The path to the desired directory in which to set up a new
#'   R package. The directory name is assumed to be the package name.
#' @return None
#' @examples
#' path <- file.path(tempdir(), "myPkgName")
#' use_ompBAM(path)
#'
#' @md
#' @export
use_ompBAM <- function(path = ".") {
    # The following packages are only required for developers using ompBAM API
    # and are added as 'suggests' to avoid having as requirements on the
    # created packages
    .check_package_is_installed("usethis")
    .check_package_is_installed("desc")
    
    proj_path <- .check_ompBAM_path(path)
    proj_name <- basename(proj_path)   
    
    usethis::create_package(proj_path, open = FALSE)
    usethis::proj_set(proj_path, force = TRUE)
    
    makevars <- file.path(proj_path, "src", "Makevars.in")
    makevars_win <- file.path(proj_path, "src", "Makevars.win")
    configure <- file.path(proj_path, "configure")
    configure.win <- file.path(proj_path, "configure.win")
    example_code <- file.path(proj_path, "src/ompBAM_example.cpp")
    R_import_file <- file.path(proj_path, "R/ompBAM_imports.R")
    usethis::use_directory("src")
    usethis::use_git_ignore(c("*.o", "*.so", "*.dll"), "src")

    .omp_write_to_package_R_script(R_import_file, proj_name)

    file.copy(system.file(file.path('examples', 'ompBAMExample', "src", 
        "Makevars.in"), package = 'ompBAM'), makevars)
    file.copy(system.file(file.path('examples', 'ompBAMExample', "src", 
        "Makevars.win"), package = 'ompBAM'), makevars_win)
    usethis::ui_done("Created src/Makevars.in and src/Makevars.win")
    file.copy(system.file(file.path('examples','ompBAMExample',"configure"), 
            package = 'ompBAM'), configure)
    file.copy(system.file(file.path('examples','ompBAMExample',"configure.win"),
            package = 'ompBAM'), configure.win)
    usethis::ui_done("Created configure scripts")
    file.copy(system.file(file.path('extdata', 'ompBAM_example.cpp'), 
        package = 'ompBAM'),example_code)
    usethis::ui_done(
        "Created src/ompBAM_example.cpp with idxstats_pbam() function")
    
    .omp_use_dependency("Rcpp", "Imports", proj_path)
    .omp_use_dependency("zlibbioc", "Imports", proj_path)
    .omp_use_dependency("ompBAM", "LinkingTo", proj_path)
    .omp_use_dependency("Rcpp", "LinkingTo", proj_path)
    .omp_use_dependency("zlibbioc", "LinkingTo", proj_path)
    
    end_msg <- paste(proj_name, "successfully created in", proj_path)
    usethis::ui_done(end_msg)
    end_msg2 <-paste("Please run roxygen2 using devtools::document() before",
        "building the package.")
    message(end_msg2)
}

#' ompBAM Example
#'
#' Installs the ompBAMExample package included with the ompBAM package
#'
#' @details
#' This function installs the ompBAMExample package located inside
#' the 'examples' subfolder of ompBAM. It uses devtools::load_all() to simulate
#' package creation. After running this function to compile ompBAMExample(), 
#' users can run the example functions contained therein.
#' @return None.
#' @examples
#' # The directory containing the source code is given by the path here
#' 
#' print(system.file(file.path('examples', "ompBAMExample"), 
#'     package = 'ompBAM'))
#'
#' # Install the ompBAMExample package
#' 
#' install_ompBAM_example()
#' @export
install_ompBAM_example <- function() {
    # The following packages are only required for developers using ompBAM API
    # and are added as 'suggests' to avoid having as requirements on the
    # created packages
    .check_package_is_installed("devtools")
    
    pkg <- "ompBAMExample"
    from <- system.file(file.path('examples', pkg), package = 'ompBAM')
    if(!dir.exists(from)) {
        stop("Invalid example package name")
    }
    dir <- tempfile()
    dir.create(dir)
    file.copy(from, dir, recursive = TRUE)
    path <- file.path(dir, pkg)

    devtools::load_all(path)
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
    dataset <- match.arg(dataset)
    if(dataset == "Unsorted") return(
        system.file(file.path('extdata', 'THP1_ND_1.bam'), package = 'ompBAM'))
    if(dataset == "scRNAseq") return(
        system.file(file.path('extdata', 'MGE_E35_1.bam'), package = 'ompBAM'))
    stop("dataset must be one of 'Unsorted' or 'scRNAseq'")
}

# Sanity checks on provided path for smooth package creation.
.check_ompBAM_path <- function(path) {   
    if(!dir.exists(dirname(path))) {
        errormsg <- paste(dirname(path), "needs to exist")
        stop(errormsg, call. = FALSE)
    }
    if(dir.exists(path)) {
        proj_path <- normalizePath(path)
    } else {
        proj_path <- file.path(normalizePath(dirname(path)), basename(path))
    }
    makevars <- file.path(proj_path, "src", "Makevars")
    makevars_win <- file.path(proj_path, "src", "Makevars.win")
    if(file.exists(makevars) | file.exists(makevars_win)) {
        fail_make_msg <- paste(
            "src/Makevars or src/Makevars.win already exist.",
            "use_ompBAM not run()"
        )
        stop(fail_make_msg, call. = FALSE)
    }
    example_code <- file.path(proj_path, "src/ompBAM_example.cpp")
    if(file.exists(example_code)) {
        fail_src_msg <- paste(
            example_code, "already exist.",
            "use_ompBAM not run()"
        )
        stop(fail_src_msg, call. = FALSE)
    }
    R_import_file <- file.path(proj_path, "R/ompBAM_imports.R")
    if(file.exists(R_import_file)) {
        fail_src_msg <- paste(
            R_import_file, "already exist.",
            "Has use_ompBAM() already been run?"
        )
        stop(fail_src_msg, call. = FALSE)
    }
    return(proj_path)
}

.omp_write_to_package_R_script <- function(R_import_file, proj_name) {
    write(paste0("#' @useDynLib ", proj_name, ", .registration = TRUE"),
        R_import_file, sep="\n")
    write("#' @import Rcpp", R_import_file, append = TRUE, sep="\n")
    write("#' @import zlibbioc", R_import_file, append = TRUE, sep="\n")
    write("NULL", R_import_file, append = TRUE, sep="\n")
    write("", R_import_file, append = TRUE, sep="\n")
    write("#' @export", R_import_file, append = TRUE, sep="\n")
    write("idxstats <- function(bam_file, n_threads) {", 
        R_import_file, append = TRUE, sep="\n")
    write("    idxstats_pbam(bam_file, n_threads)", 
        R_import_file, append = TRUE, sep="\n")
    write("}", R_import_file, append = TRUE, sep="\n")

    usethis::ui_done(paste("Created", R_import_file, "with roxygen tags"))
}

# modeled after usethis::use_dependency
.omp_use_dependency <- function(package, type, proj_path) {
    types <- c("Depends", "Imports", "Suggests", "Enhances", "LinkingTo")
    names(types) <- tolower(types)
    type <- types[[match.arg(tolower(type), names(types))]]
    
    deps <- desc::desc_get_deps(proj_path)
    
    existing_dep <- deps$package == package
    existing_type <- deps$type[existing_dep]
    is_linking_to <- (existing_type != "LinkingTo" & type == "LinkingTo") |
        (existing_type == "LinkingTo" & type != "LinkingTo")
    
    if (!any(existing_dep) || any(is_linking_to)) {
        done_msg <- paste("Adding", package, "to", type, "field in DESCRIPTION")
        usethis::ui_done(done_msg)
        desc::desc_set_dep(package, type, file = proj_path)
        return(invisible(TRUE))
    }
    
    existing_type <- setdiff(existing_type, "LinkingTo")
    delta <- sign(match(existing_type, types) - match(type, types))
    if (delta <= 0) {
        # don't downgrade
        warn_msg <- paste(
            "Package", package, "is already listed in",
            existing_type, "in DESCRIPTION, no change made."
        )
        usethis::ui_warn(warn_msg)
        return(invisible(FALSE))
    } else if (delta > 0) {
    # upgrade
        if (existing_type != "LinkingTo") {
            done_msg <- paste("Moving", package, "from", existing_type,
                "to", type, "field in DESCRIPTION")
            usethis::ui_done(done_msg)
            desc::desc_del_dep(package, existing_type, file = proj_path)
            desc::desc_set_dep(package, type, file = proj_path)
        }
    }
    invisible(TRUE)
}

.check_package_is_installed <- function(
        package = "DESeq2", 
        version = "1.0.0", 
        returntype = c("error", "warning", "silent")
) {
    res <- tryCatch(
        ifelse(utils::packageVersion(package)>=version, TRUE, FALSE),
        error = function(e) FALSE)
    if(!res) {
        returntype <- match.arg(returntype)
        stopmsg <- paste(package, "version", version, "is not installed;",
            "and is required for this function")
        if(returntype == "error") {
            stop(stopmsg, call. = FALSE)
        } else if(returntype == "warning") {
            warning(stopmsg, call. = FALSE)
        } else {
            message(stopmsg)
        }
    }
    return(res)
}