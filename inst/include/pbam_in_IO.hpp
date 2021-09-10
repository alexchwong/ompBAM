#ifndef _pbam_in_IO
#define _pbam_in_IO

// Public functions:

inline int pbam_in::openFile(std::string filename, unsigned int n_threads) {
  check_threads(n_threads);
  clear_buffers();
  
  IN = new std::ifstream(filename, std::ios::in | std::ifstream::binary);
  FILENAME = filename;
  int ret = check_file();
  return(ret);
}

inline int pbam_in::SetInputHandle(std::ifstream *in_stream, unsigned int n_threads) {
  check_threads(n_threads);
  if(!in_stream) return(-1);
  
  IN = in_stream;
  int ret = check_file();
  return(ret);
}

inline int pbam_in::closeFile() {
  clear_buffers(); return(0);
}

inline int pbam_in::obtainChrs(std::vector<std::string> & s_chr_names, std::vector<uint32_t> & u32_chr_lens) {
  if(!magic_header) {
    Rcpp::Rcout << "Header is not yet read\n";
    return(-1);
  }
  if(n_ref == 0) {
    Rcpp::Rcout << 
      "No chromosome names stored. Is pbam_in::readHeader() been run yet?\n";
    return(-1);
  }
  s_chr_names.clear();
  u32_chr_lens.clear();
  for(unsigned int i = 0; i < n_ref; i++) {
    s_chr_names.push_back(chr_names.at(i));
    u32_chr_lens.push_back(chr_lens.at(i));
  }
  return((int)n_ref);
}

// Internals

inline int pbam_in::check_file() {
  if(!IN->fail()) {
    // Assign IS_LENGTH
    IN->seekg(0, std::ios_base::end);
    IS_LENGTH = tellg();
    
    // Check valid BAM EOF:
    IN->seekg(IS_LENGTH-bamEOFlength, std::ios_base::beg);
    
    char eof_check[bamEOFlength + 1];
    IN->read(eof_check, bamEOFlength);
    if(strncmp(bamEOF, eof_check, bamEOFlength) != 0) {
      Rcpp::Rcout << "Error opening BAM - EOF bit corrupt. Perhaps this file is truncated?\n";
      IN = NULL;
      return(-1);
    }
    IN->clear(); 
    IN->seekg(0, std::ios_base::beg);
    
    // Read header. If corrupt header, close everything and output error:
    int ret = readHeader();
    if(ret != 0) {
      clear_buffers();
    }
    return(ret);
     // return(0);
  } else {
    return(-1);
  }
}

inline int pbam_in::readHeader() {
  if(magic_header) {
    Rcpp::Rcout << "Header is already read\n";
    return(-1);
  }
  
  magic_header = (char*)malloc(8+1);
  read(magic_header, 8);
  if(strncmp(magic_header, magicstring, 4) != 0) {
    Rcpp::Rcout << "Invalid BAM magic string\n";
    free(magic_header); magic_header = NULL;
    return(-1);
  }
  
  uint32_t * u32 = (uint32_t *)(magic_header + 4);
  l_text = *u32;
  headertext = (char*)malloc(l_text + 1);
  read(headertext, l_text);
  
  char * u32c = (char*)malloc(5);
  read(u32c, 4);
  u32 = (uint32_t *)u32c;
  n_ref = *u32;
  char chrom_buffer[1000]; std::string chrName;
  for(unsigned int i = 0; i < n_ref; i++) {
    read(u32c, 4); u32 = (uint32_t *)u32c;
    read(chrom_buffer, *u32);
    chrName = std::string(chrom_buffer, *u32-1);
    chr_names.push_back(chrName);

    read(u32c, 4); u32 = (uint32_t *)u32c;
    chr_lens.push_back(*u32);
  }
  
  free(u32c);
  return(0);
}

#endif
