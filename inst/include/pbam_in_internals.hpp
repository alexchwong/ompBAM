#ifndef _pbam_in_internals
#define _pbam_in_internals

// All functions in this file are private functions.

inline void pbam_in::initialize_buffers() {
  // Empty buffer pointers
  file_buf = NULL;  data_buf = NULL;  next_file_buf = NULL;
  // Empty capacities and cursors
  file_buf_cap = 0; file_buf_cursor = 0;
  data_buf_cap = 0; data_buf_cursor = 0;
  next_file_buf_cap = 0; next_file_buf_cursor = 0;
  // Empty BAM header info
  magic_header = NULL; l_text = 0; headertext = NULL; n_ref = 0;
  chr_names.resize(0); chr_lens.resize(0);
  // Clears file name string
  FILENAME.clear();

  // Empties cursors for thread-specific reads
  read_cursors.resize(0); read_ptr_ends.resize(0);

  // Clears handle to ifstream
  IN = NULL;
  error_state = 0;
}

inline void pbam_in::clear_buffers() {
  if(FILENAME.size() > 0 && IN) {
    IN->close();  // close previous file if opened using openFile()
    delete(IN);   // Releases heap-allocated ifstream produced by openFile()
    FILENAME.clear();
  }
  
  // Releases file and data buffers
  if(file_buf) free(file_buf); 
  file_buf = NULL;
  if(data_buf) free(data_buf); 
  data_buf = NULL;
  if(next_file_buf) free(next_file_buf); 
  next_file_buf = NULL;
  file_buf_cap = 0; file_buf_cursor = 0;
  data_buf_cap = 0; data_buf_cursor = 0;
  next_file_buf_cap = 0; next_file_buf_cursor = 0;

  // Releases stored header data
  if(magic_header) free(magic_header); 
  magic_header = NULL;
  if(headertext) free(headertext); 
  headertext = NULL;  
  l_text = 0; n_ref = 0;
  chr_names.resize(0); chr_lens.resize(0);

  // Empties cursors for thread-specific reads
  read_cursors.resize(0);
  read_ptr_ends.resize(0);

  // Clears handle to ifstream
  IN = NULL;
}

// Makes sure given threads does not exist system resources
inline void pbam_in::check_threads(unsigned int n_threads_to_check) {
  #ifdef _OPENMP
    if(n_threads_to_check > (unsigned int)omp_get_max_threads()) {
      threads_to_use = (unsigned int)omp_get_max_threads();
    } else {
      threads_to_use = n_threads_to_check;
    }
  #else
    threads_to_use = 1;
  #endif
}

// Analogous functions to fstream, but for reading from decompressed buffer

// Reads bytes from the current data_buf, increments data_buff_cursor
inline unsigned int pbam_in::read(char * dest, const unsigned int len) {
  if(data_buf_cap - data_buf_cursor < len) {
    decompress(len + 65536);
  }
  
  unsigned int n_bytes_to_read = std::min((size_t)len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}

// Ignores bytes from the current data_buf, increments data_buff_cursor
inline unsigned int pbam_in::ignore(const unsigned int len) {
  if(data_buf_cap - data_buf_cursor < len) {
    decompress(len + 65536);
  }
  
  unsigned int n_bytes_to_read = std::min((size_t)len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  // memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}

// Peeks at bytes from the current data_buf, does not increment data_buff_cursor
inline unsigned int pbam_in::peek(char * dest, const unsigned int len) {
  if(data_buf_cap - data_buf_cursor < len) {
    decompress(len + 65536);
  }
  
  unsigned int n_bytes_to_read = std::min((size_t)len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  // data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}


#endif
