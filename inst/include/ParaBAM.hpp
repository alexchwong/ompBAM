#ifndef _ParaBAM
#define _ParaBAM

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
using namespace Rcpp;

#include "pbam_defs.hpp"

#include "pbam1_t.hpp"
#include "pbam_in.hpp"

#endif

