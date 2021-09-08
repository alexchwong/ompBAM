// [[Rcpp::depends(ompBAM)]]
#include <ompBAM.hpp>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
int idxstats_pbam(std::string bam_file, int n_threads_to_use = 1, bool verbose = true){

  unsigned int n_threads_to_really_use;
  #ifdef _OPENMP
    if(n_threads_to_use > 1) {
      if(n_threads_to_use > omp_get_max_threads()) {
        n_threads_to_really_use = omp_get_max_threads();
      } else {
        n_threads_to_really_use = n_threads_to_use;
      }
    } else {
      n_threads_to_really_use = 1;
    }
  #else
    n_threads_to_really_use = 1;
  #endif

  pbam_in inbam;
  inbam.openFile(bam_file, n_threads_to_really_use);
  
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  int chrom_count = inbam.obtainChrs(s_chr_names, u32_chr_lens);

  if(chrom_count <= 0) {
    return(-1); // obtainChrs and SetInputHandle already returns relevant error msgs
  }
  
  // Creates a data structure that stores per-chromosome read counts
  std::vector<uint32_t> total_reads(chrom_count);

  Progress p(inbam.GetFileSize(), verbose);
  bool report_first_read = false;
  while(0 == inbam.fillReads()) {
    // Rcout << "Reads filled\n";
    p.increment(inbam.IncProgress());
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_really_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_really_use; i++) {
      // Rcout << "Thread " << i << "\n";
      std::vector<uint32_t> read_counter(chrom_count);
      pbam1_t read;
      read = inbam.supplyRead(i);
      do {
        if(read.validate()) {
          if(read.refID() >= 0) {
            read_counter.at(read.refID())++;
            
            // Check read getters
            if(!report_first_read && i == 0) {
              std::string read_name;
              read.read_name(read_name);
              Rcout << "Read name: " << read_name << '\n';
              
              std::string read_seq;
              read.seq(read_seq);
              Rcout << "Seq: " << read_seq << '\n';

              std::string cigar;
              read.cigar(cigar);
              Rcout << "Cigar: " << cigar << '\n';
              
              Rcout << "Available tags: ";
              std::vector<std::string> tags;
              read.AvailTags(tags);
              for(unsigned int j = 0; j < tags.size(); j++) {
                Rcout << tags.at(j) << " ";
              }
              
              Rcout << '\n';
              #ifdef _OPENMP
              #pragma omp critical
              #endif
              report_first_read = true;
            }            
          }
        }
        read = inbam.supplyRead(i);
      } while(read.validate());
    
      // Summarise reads:
      #ifdef _OPENMP
      #pragma omp critical
      #endif
      for(unsigned int j = 0; j < (unsigned int)chrom_count; j++) {
        total_reads.at(j) += read_counter.at(j);
      }
    }
  }
  p.increment(inbam.IncProgress());
  
  // inbam_stream.close();

  Rcout << bam_file << " summary:\n" << "Name\tLength\tNumber of reads\n";
  for(unsigned int j = 0; j < (unsigned int)chrom_count; j++) {
    Rcout << s_chr_names.at(j) << '\t' << u32_chr_lens.at(j) << '\t'
      << total_reads.at(j) << '\n';
  }
  return(0);
}