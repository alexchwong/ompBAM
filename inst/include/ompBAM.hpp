#ifndef _ompBAM
#define _ompBAM

#include <fstream>    // std::ifstream

// [[Rcpp::depends(zlibbioc)]]
#include <zlib.h>     // For BAM decompression
#include <zconf.h>

#include <cstring>    // To compare between strings
#include <vector>     // For vector types

#ifdef _OPENMP
  #include <omp.h>    // For OpenMP
#endif

#include "Rcpp.h"     // For Rcpp::Rcout

#include "pbam_defs.hpp"
#include "pbam1_t.hpp"
#include "pbam_in.hpp"

inline void paraBAM_version() {
  std::string version = "0.1.0";
  Rcpp::Rcout << "Compiled using ompBAM version " << version << "\n";
}

#endif

