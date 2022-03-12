/* ompBAM.hpp ompBAM Library

Copyright (C) 2021 Alex Chit Hei Wong

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.  */

#ifndef _ompBAM
#define _ompBAM

#include <fstream>    // std::ifstream

// [[Rcpp::depends(zlibbioc)]]
#include <zlib.h>     // For BAM decompression
#include <zconf.h>

#include <string>    
#include <cstring>    // To compare between strings
#include <vector>     // For vector types
#include <iostream>   // For cout

#ifdef _OPENMP
  #include <omp.h>    // For OpenMP
#endif

#include "pbam_defs.hpp"
#include "pbam1_t.hpp"
#include "pbam_in.hpp"

inline void ompBAM_version() {
  std::string version = "0.99.0";
  cout << "Compiled using ompBAM version " << version 
    << " for NxtIRFcore\n";
}

inline void is_compiled_with_OpenMP() {
  #ifdef _OPENMP
    cout << "Compiled with OpenMP; " << omp_get_max_threads()
      << " threads available\n";
  #else
    cout << "Compiled without OpenMP support\n";
  #endif
}

#endif

