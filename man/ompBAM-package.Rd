% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ompBAM.R
\docType{package}
\name{ompBAM-package}
\alias{ompBAM}
\alias{ompBAM-package}
\title{ompBAM: C++ header library for parallel sequential reading
of BAM files}
\description{
Although various tools for multi-threaded BAM processing (samtools,
sambamba) are available, currently there are none available to the
R/Bioconductor environment. Parallelism in R is achieved using BiocParallel
or related packages. Typical alternatives include processing multiple BAM
files, each using a single core. Although easy to set up, memory-intensive
applications limit such approaches. ompBAM is a header-only C++ library
based on OpenMP, allowing developers to implement OpenMP-based parallelism
for the sequential reading of BAM files. ompBAM makes it easy by handling
all the parallel file-access and bgzf decompression, allowing developers to
focus on multi-threaded handling of individual reads.
}
\details{
For a quick start guide and a full list of functionality, refer to
the vignette
\href{../doc/ompBAM-API-Docs.html}{\code{vignette("ompBAM-API-Docs", 
package = "ompBAM")}}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/alexchwong/ompBAM}
  \item Report bugs at \url{https://support.bioconductor.org/}
}

}
\author{
Alex Chit Hei Wong
}
\keyword{package}
