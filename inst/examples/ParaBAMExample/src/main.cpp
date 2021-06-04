// [[Rcpp::depends(ParaBAM)]]
#include <ParaBAM.hpp>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
int idxstats_pbam(std::string bam_file, int n_threads_to_use = 1, bool verbose = true){

  unsigned int n_threads_to_really_use;
  #ifdef _OPENMP
    n_threads_to_really_use = std::max(n_threads_to_use, 1);
  #else
    n_threads_to_really_use = 1;
  #endif

  
  std::ifstream inbam_stream;   
  inbam_stream.open(bam_file, std::ios::in | std::ios::binary);

  // File buffer 1 Gb, Data buffer 2 Gb, 10 file chunks per buffer
  // Buffer usage = 2 * 1 Gb (File + nextFile buffers) + 2Gb (Data buffer)
  pbam_in inbam((size_t)1000000000, (size_t)2000000000, 10);
  inbam.SetInputHandle(&inbam_stream, n_threads_to_really_use);
  
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  int ret = inbam.obtainChrs(s_chr_names, u32_chr_lens);

  if(ret <= 0) {
    return(-1); // obtainChrs and SetInputHandle already returns relevant error msgs
  }
  
  // Creates a data structure that stores per-chromosome read counts
  std::vector<uint32_t> total_reads(ret);

  Progress p(inbam.GetFileSize(), verbose);

  while(0 == inbam.fillReads()) {
    p.increment(inbam.IncProgress());
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_really_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_really_use; i++) {
      std::vector<uint32_t> read_counter(ret);
      pbam1_t read;
      read = inbam.supplyRead(i);
      do {
        if(read.validate()) {
          if(read.refID() >= 0) {
            read_counter.at(read.refID())++;
          }
        }
        read = inbam.supplyRead(i);
      } while(read.validate());
    
      // Summarise reads:
      #ifdef _OPENMP
      #pragma omp critical
      #endif
      for(unsigned int j = 0; j < (unsigned int)ret; j++) {
        total_reads.at(j) += read_counter.at(j);
      }
    }
  }

  inbam_stream.close();

  Rcout << bam_file << " summary:\n" << "Name\tLength\tNumber of reads\n";
  for(unsigned int j = 0; j < (unsigned int)ret; j++) {
    Rcout << s_chr_names.at(j) << '\t' << u32_chr_lens.at(j) << '\t'
      << total_reads.at(j) << '\n';
  }
  return(0);
}