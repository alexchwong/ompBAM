#ifndef _ompBAM
#define _ompBAM

#include <fstream>    // std::ifstream

// [[Rcpp::depends(zlibbioc)]]
#include <zlib.h>
#include <zconf.h>

#include <cstring>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "Rcpp.h"

#include "pbam_defs.hpp"

#include "pbam1_t.hpp"
#include "pbam_in.hpp"

inline void paraBAM_version() {
  std::string version = "0.1.0";
  Rcpp::Rcout << "Compiled using ompBAM version " << version << "\n";
}

#endif

