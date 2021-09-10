# ompBAM
C++ Library for OpenMP-based multi-threaded sequential profiling of Binary Alignment Map (BAM) files

# Description
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

## Installation

### On current R (>= 4.0.0)
* Development version from Github:
```
library("devtools")
install_github("alexchwong/ompBAM", dependencies=TRUE)
```

## Documentation

Access the ompBAM-API-Docs via its included vignette. This includes:
* How to set up a new package R-project, ready-to-compile with ompBAM, as well as a 'Hello World' equivalent example function of the 'idxstats' function to demonstrate ompBAM
* A step-by-step guide of how the idxstats function implemented in the example code is constructed
* Detailed documentation of the `pbam_in` and `pbam1_t` objects that comprise ompBAM.

```
browseVignettes("NxtIRF")
```
