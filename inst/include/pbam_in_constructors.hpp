#ifndef _pbam_in_constructors
#define _pbam_in_constructors

// Declare pbam_in with default settings
inline pbam_in::pbam_in() {
  initialize_buffers();
}

// Declare pbam_in with custom settings; for advanced users only
inline pbam_in::pbam_in(
  const size_t file_buffer_cap,   // Default 500 Mb
  const size_t data_buffer_cap,   // Default 1 Gb
  const unsigned int chunks_per_file_buffer,    // Default 5 - i.e. 40 Mb file chunks
  // File read triggers after each 40 Mb decompressed
  const bool read_file_using_multiple_threads   // default true
) {
  initialize_buffers();
  
  // Limits chunk size to be above 1 Mb
  if(file_buffer_cap / chunks_per_file_buffer < 1024576) {
    Rcpp::Rcout << "FILE_BUFFER_CAP / chunks_per_file_buffer (chunk size) must be above 1Mb\n";
    return;
  }
  // Ensures decompressed data buffer is bigger than file buffer. For sanity reasons.
  if(data_buffer_cap < file_buffer_cap) {
    Rcpp::Rcout << "DATA_BUFFER_CAP must not be smaller than FILE_BUFFER_CAP\n";
    return;
  }  

  FILE_BUFFER_CAP = file_buffer_cap;
  DATA_BUFFER_CAP = data_buffer_cap;
  chunks_per_file_buf = chunks_per_file_buffer;
  threads_to_use = 1;
  multiFileRead = read_file_using_multiple_threads;
}

inline pbam_in::~pbam_in() {
  clear_buffers();
}

#endif
