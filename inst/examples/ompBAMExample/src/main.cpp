// [[Rcpp::depends(ompBAM)]]
#include <ompBAM.hpp>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

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
      pbam1_t read(inbam.supplyRead(i));
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
              
              std::vector<std::string> tags;
              read.AvailTags(tags);
              for(unsigned int j = 0; j < tags.size(); j++) {
                std::string Z_val;
                Rcout << tags.at(j) << ":" << read.Tag_Type_SAM(tags.at(j)) << ":";
                switch(read.Tag_Type(tags.at(j))) {
                  case 'A':
                    Rcout << read.tagVal_A(tags.at(j)) << '\t'; break; 
                  case 'c':
                    Rcout << std::to_string(read.tagVal_c(tags.at(j))) << '\t'; break; 
                  case 'C':
                    Rcout << std::to_string(read.tagVal_C(tags.at(j))) << '\t'; break; 
                  case 's':
                    Rcout << std::to_string(read.tagVal_s(tags.at(j))) << '\t'; break; 
                  case 'S':
                    Rcout << std::to_string(read.tagVal_S(tags.at(j))) << '\t'; break; 
                  case 'i':
                    Rcout << std::to_string(read.tagVal_i(tags.at(j))) << '\t'; break; 
                  case 'I':
                    Rcout << std::to_string(read.tagVal_I(tags.at(j))) << '\t'; break; 
                  case 'f':
                    Rcout << std::to_string(read.tagVal_f(tags.at(j))) << '\t'; break; 
                  case 'Z':
                    read.tagVal_Z(tags.at(j), Z_val);
                    Rcout << Z_val << '\t'; break;
                }
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