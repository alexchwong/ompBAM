#include <Rcpp.h>
#include "ParaBAM.h"

// [[Rcpp::export]]
int main(std::string bam_file){

  std::ifstream inbam_stream;   
  inbam_stream.open(bam_file, std::ios::in | std::ios::binary);

  ParaBAM inbam;
  inbam.SetInputHandle(&inbam_stream);
  
  inbam.readHeader();
  
  std::vector<std::string> & s_chr_names;
  std::vector<uint32_t> & u32_chr_lens;
  
  int ret = inbam.obtainChrs(s_chr_names, u32_chr_lens);
  
  for(unsigned int i = 0; i < ret; i++) {
    Rcout << s_chr_names.at(i) << '\t' << u32_chr_lens.at(i) << '\n';
  }
  return(0);
}