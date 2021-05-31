#include "ParaBAM.h"
// [[Rcpp::depends(ParaBAM)]]

#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
int test_ParaBAM(std::string bam_file, int n_threads_to_use = 1){

  std::ifstream inbam_stream;   
  inbam_stream.open(bam_file, std::ios::in | std::ios::binary);

  ParaBAM inbam;
  inbam.SetInputHandle(&inbam_stream, (unsigned int)n_threads_to_use);
  
  inbam.readHeader();
  
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  
  int ret = inbam.obtainChrs(s_chr_names, u32_chr_lens);
  
  for(int i = 0; i < ret; i++) {
    Rcout << s_chr_names.at(i) << '\t' << u32_chr_lens.at(i) << '\n';
  }
  
  inbam.fillReads();
  
  uint32_t total_reads = 0;
  
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(unsigned int i = 0; i < (unsigned int)n_threads_to_use; i++) {
    unsigned int read_counter = 0;
    char * read;
    do {
      read = inbam.supplyRead(i);
      if(read) {
        read_counter++;

        // Prints the read names of first 10 reads of each thread
        if(read_counter <= 10) {
          uint8_t l_read_name;
          char * read_name = inbam.readName(read, l_read_name);
          std::string s_read_name(read_name, l_read_name);
          
          #ifdef _OPENMP
          #pragma omp critical
          #endif
          Rcout << "thread " << i << " read " << s_read_name << '\n';
        }

      }
    } while(read);
  
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    total_reads += read_counter;
  }
  inbam_stream.close();

  Rcout << "Total reads = " << total_reads << '\n';
  return(0);
}