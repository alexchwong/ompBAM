#include <ompBAM.hpp>

#include "Rcpp.h"
using namespace Rcpp;

unsigned int use_threads(int n_threads = 1) {
  #ifdef _OPENMP
  if(n_threads <= 1) return(1);
  if(n_threads > omp_get_max_threads()) return((unsigned int)omp_get_max_threads());
  return((unsigned int)n_threads);
  #else
  return(1);
  #endif
}

// [[Rcpp::export]]
int idxstats_pbam(std::string bam_file, int n_threads_to_use = 1, bool verbose = true){
  
  // Ensure number of threads requested < number of system threads available
  unsigned int n_threads_to_really_use = use_threads(n_threads_to_use);
  Rcout << "Using " << n_threads_to_really_use << " threads\n";

  // Open BAM file
  pbam_in inbam;
  inbam.openFile(bam_file, n_threads_to_really_use);
  
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  int chrom_count = inbam.obtainChrs(s_chr_names, u32_chr_lens);

  // obtainChrs and SetInputHandle already returns relevant error msgs
  if(chrom_count <= 0) return(-1);
  
  // Creates a data structure that stores per-chromosome read counts
  std::vector<uint32_t> total_reads(chrom_count);

  // Loop to read part BAM file, filling the buffers with reads
  while(0 == inbam.fillReads()) {
    // Internal OpenMP parallel FOR loop, 1 thread runs 1 loop, n threads
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_really_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_really_use; i++) {
      std::vector<uint32_t> read_counter(chrom_count);
      pbam1_t read(inbam.supplyRead(i));
      do {
        if(read.validate()) {
          if(read.refID() >= 0) {
            read_counter.at(read.refID())++;     
          }
        }
        read = inbam.supplyRead(i);
      } while(read.validate());
    
      // Add reads tally. pragma omp critical ensures only 1 thread writes at a time
      #ifdef _OPENMP
      #pragma omp critical
      #endif
      for(unsigned int j = 0; j < (unsigned int)chrom_count; j++) {
        total_reads.at(j) += read_counter.at(j);
      }
    }
  }

  // Outputs summarised reads
  Rcout << bam_file << " summary:\n" << "Name\tLength\tNumber of reads\n";
  for(unsigned int j = 0; j < (unsigned int)chrom_count; j++) {
    Rcout << s_chr_names.at(j) << '\t' << u32_chr_lens.at(j) << '\t'
      << total_reads.at(j) << '\n';
  }
  return(0);
}