# From rcpp_progress

load_my_example_pkg <- function(pkg, ...) {
  skip_if(!requireNamespace("devtools", quietly = TRUE),
    message = "Package devtools must be installed to run unit tests.")

  from <- system.file(file.path('examples', pkg), package = 'ParaBAM')
  dir <- tempfile()
  dir.create(dir)
  file.copy(from, dir, recursive = TRUE)
  path <- file.path(dir, pkg)

  devtools::load_all(path, quiet = TRUE, ...)
}

get_function_from_pkg <- function(pkg, fun) {
  get(fun, getNamespace(pkg))
}


test_example <- function(...) {
  pkg <- 'ParaBAMExample'
  load_my_example_pkg(pkg, ...)
  fun <- get_function_from_pkg(pkg, 'idxstats_pbam')
  
  bam <- system.file(file.path('extdata', 'THP1.bam'), package = 'ParaBAM')
  
  fun(bam, n_threads_to_use = 2)
}