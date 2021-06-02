
install_ParaBAM_example <- function(pkg = "ParaBAMExample") {
    from <- system.file(file.path('examples', pkg), package = 'ParaBAM')
    if(!dir.exists(from)) {
        stop("Invalid example package name")
    }
    dir <- tempfile()
    dir.create(dir)
    file.copy(from, dir, recursive = TRUE)
    path <- file.path(dir, pkg)

    devtools::load_all(path, quiet = TRUE, ...)
}