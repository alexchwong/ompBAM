#ifndef _pbam_in
#define _pbam_in

/*
  Class Description
*/
class pbam_in {
  public:
    // Creates pbam_in using default values
    pbam_in();
    
    // Creates a pbam_in with custom caps on file and data buffer sizes, 
    // as well as chunks per file
    pbam_in(
      const size_t file_buffer_cap, 
      const size_t data_buffer_cap, 
      const unsigned int chunks_per_file_buffer,
      const bool read_file_using_multiple_threads = true
    );
    
    ~pbam_in();

    /*
      Opens a file directly for BAM reading. Supply a file name as a string
    */
    int openFile(std::string filename, unsigned int n_threads); 
    int closeFile();
    
    /*
      Assign a ifstream handle to pbam_in; declares number of threads to use.
      in_stream must point to an ifstream that has opened a BAM file in binary mode
      Additionally, this function reads the BAM header.
    */
    int SetInputHandle(std::ifstream *in_stream, const unsigned int n_threads); 
    
    // Returns two vectors, containing chromosome names and lengths
    int obtainChrs(
      std::vector<std::string> & s_chr_names, 
      std::vector<uint32_t> & u32_chr_lens
    );

    /* 
      Reads the BAM file, decompressing to a maximum either by the data buffer cap,
        or the file chunk size (file_buf_cap / chunks_per_file_buf).
      Then assigns read cursors to point to the locations of reads to begin reading
        by each thread.
    */
    int fillReads();  // Returns 1 if no more reads available
    
    /* 
      Returns a pbam1_t of the next read to be returned by thread 
        (as specified by thread_number).
      - thread_number must be between [0, n_threads - 1].
      - If there are no more reads in the thread, returns an empty read which will
        return false using pbam1_t::validate()
    */
    pbam1_t supplyRead(const unsigned int thread_number = 0);
    
    // Returns the size of the opened BAM
    size_t GetFileSize() { return(IS_LENGTH); };

    // Returns the number of bytes decompressed
    size_t GetProgress() {return(prog_tellg());};
    
    /* 
      Returns the incremental number of bytes decompressed since the last call 
      to IncProgress()
    */
    size_t IncProgress() {
      size_t INC = prog_tellg() - PROGRESS;
      PROGRESS = prog_tellg();
      return(INC);
    };

  private:
// pbam_in Settings
    size_t          FILE_BUFFER_CAP       = 2e8;
    size_t          DATA_BUFFER_CAP       = 1e9;
    unsigned int    chunks_per_file_buf   = 5;    // Divide file buffer into n segments
    unsigned int    threads_to_use        = 1;
    bool            multiFileRead         = true;
    std::string     FILENAME;
// File particulars
    std::ifstream    * IN;    
    size_t          IS_LENGTH;     // size of BAM file (Input Stream Length)
  
// Header storage:
    char                        * magic_header;  // Always of size = 8
    uint32_t                    l_text;      // sizeof *headertext
    char                        * headertext;
    
// Chromosomes:
    uint32_t                    n_ref;    
    std::vector<std::string>    chr_names;
    std::vector<uint32_t>       chr_lens;
  
// Data buffers
    char *          file_buf; 
    size_t          file_buf_cap; 
    size_t          file_buf_cursor;
    
    char *          next_file_buf; 
    size_t          next_file_buf_cap;
    
    char *          data_buf; 
    size_t          data_buf_cap; 
    size_t          data_buf_cursor;
    
/* 
  Thread-specific read cursor positions and boundaries
  Thread returns a null read if cursor read_cursors >= read_ptr_ends
*/
    std::vector<size_t>         read_cursors;     // Cursor(s) of start of next read in each thread
    std::vector<size_t>         read_ptr_ends;    // Boundaries of read positions in each thread

// Internal functions

    void            check_threads(unsigned int n_threads_to_check);

    int             check_file();

// *** Initialisers ***
    void            initialize_buffers();         // Initialises a pbam_in
    void            clear_buffers();              // Clears all buffers and re-initialises pbam_in

// *** Internal functions run by decompress() ***

    int             read_file_to_buffer(char * buf, const size_t len);

    // Removes all data upstream of data_buf_cursor
    int             clean_data_buffer(const size_t n_bytes_to_decompress);  
    
    /* 
      If bytes remaining in primary file buffer is less than chunk_size,
      re-initialises file_buf, copying residual data from existing file_buf,
      and copies data from next_file_buf, upto FILE_BUFFER_CAP
    */
    int             swap_file_buffer_if_needed();
    
    // Reads from file and copies data to file_buf, upto a maximum of n_bytes
    size_t          load_from_file(const size_t n_bytes);

    // Runs load_from_file using n_bytes = FILE_BUFFER_CAP
    size_t          fill_file_buffer();

    // Essentially same as load_from_file, but this reads data into next_file_buf
    size_t          read_file_chunk_to_spare_buffer(const size_t n_bytes);

// *** Main DECOMPRESS function ***
    size_t          decompress(const size_t n_bytes_to_decompress);

// *** Internal functions used by readHeader() ***
    unsigned int read(char * dest, const unsigned int len);  // returns the number of bytes actually read
    unsigned int ignore(unsigned const int len);
    unsigned int peek(char * dest, const unsigned int len) ;
    
// *** Reads the BAM header. Automatically run with SetInputHandle() ***
    int             readHeader();


// *** File specific functions ***
    size_t tellg() {return((size_t)IN->tellg());};    // Returns position of file cursor
    
    size_t prog_tellg() {                             // Returns the number of bytes decompressed
      return((size_t)IN->tellg()
        - (file_buf_cap-file_buf_cursor)
        - next_file_buf_cap
    ); };
    
    bool eof() {return(IS_LENGTH <= tellg());};       // Returns whether end of file is reached
    bool fail() {return(IN->fail());};                // Returns any ifstream errors
    
    size_t PROGRESS = 0;    // Value of prog_tellg() when IncProgress() is last called
    
// Disable copy construction / assignment (doing so triggers compile errors)
    pbam_in(const pbam_in &t);
    pbam_in & operator = (const pbam_in &t);

};

// Functions

inline void pbam_in::initialize_buffers() {
  file_buf = NULL;
  data_buf = NULL;
  next_file_buf = NULL;
  file_buf_cap = 0; file_buf_cursor = 0;
  data_buf_cap = 0; data_buf_cursor = 0;
  next_file_buf_cap = 0;
  magic_header = NULL; l_text = 0; headertext = NULL; n_ref = 0;
  chr_names.resize(0); chr_lens.resize(0);

  read_cursors.resize(0);
  read_ptr_ends.resize(0);
  IN = NULL;
}

inline void pbam_in::clear_buffers() {
  if(FILENAME.size() > 0 && IN) {
    IN->close(); // close previous file
    delete(IN);
    FILENAME.clear();
  }
  if(file_buf) free(file_buf); 
  file_buf = NULL;
  if(data_buf) free(data_buf); 
  data_buf = NULL;
  if(next_file_buf) free(next_file_buf); 
  next_file_buf = NULL;

  if(magic_header) free(magic_header); 
  magic_header = NULL;
  if(headertext) free(headertext); 
  headertext = NULL;
  
  file_buf_cap = 0; file_buf_cursor = 0;
  data_buf_cap = 0; data_buf_cursor = 0;
  next_file_buf_cap = 0;
  l_text = 0; n_ref = 0;
  chr_names.resize(0); chr_lens.resize(0);

  read_cursors.resize(0);
  read_ptr_ends.resize(0);
  IN = NULL;
}

inline pbam_in::pbam_in() {
  initialize_buffers();
}

inline pbam_in::pbam_in(
  const size_t file_buffer_cap, 
  const size_t data_buffer_cap, 
  const unsigned int chunks_per_file_buffer,
  const bool read_file_using_multiple_threads
) {
  initialize_buffers();
  
  if(file_buffer_cap / chunks_per_file_buffer < 1024576) {
    Rcpp::Rcout << "FILE_BUFFER_CAP / chunks_per_file_buffer (chunk size) must be above 1Mb\n";
    return;
  }
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

inline int pbam_in::check_file() {
  if(IN) {
    // Assign IS_LENGTH
    IN->seekg(0, std::ios_base::end);
    IS_LENGTH = tellg();
    
    // Check valid EOF:
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
  } else {
    return(-1);
  }
}

inline int pbam_in::openFile(std::string filename, unsigned int n_threads) {
  check_threads(n_threads);
  clear_buffers();
  
  IN = new std::ifstream;
  FILENAME = filename;
  IN->open(filename, std::ios::in | std::ifstream::binary);
  return(check_file());
}

inline int pbam_in::SetInputHandle(std::ifstream *in_stream, unsigned int n_threads) {
  check_threads(n_threads);
  
  IN = in_stream;
  return(check_file());
}

inline int pbam_in::closeFile() {
  clear_buffers(); return(0);
}

inline int pbam_in::swap_file_buffer_if_needed() {
  // Transfers residual data from current file buffer to the next:
  if(next_file_buf_cap == 0) return(1);
  size_t chunk_size = (size_t)(FILE_BUFFER_CAP / chunks_per_file_buf);
  if(file_buf_cap  - file_buf_cursor > chunk_size) return(1);
  
  char * file_tmp;
  char * residual_data_buffer;
  size_t residual = file_buf_cap-file_buf_cursor;
  // Remove residual bytes to beginning of buffer:
  if(residual > 0) {
    // have to move bytes around
    residual_data_buffer = (char*)malloc(residual + 1);
    memcpy(residual_data_buffer, file_buf + file_buf_cursor, residual);
  
    file_buf = (char*)realloc(file_tmp = file_buf, FILE_BUFFER_CAP + 1);
    memcpy(file_buf, residual_data_buffer, residual);
    free(residual_data_buffer);
  } else {
    // residual == 0
    file_buf = (char*)realloc(file_tmp = file_buf, FILE_BUFFER_CAP + 1);
  }
  
  file_buf_cap = residual;
  file_buf_cursor = 0;
  
  // If next_file_buf fits inside cap, transfer all contents over and free next_file_buf
  if(next_file_buf_cap <= FILE_BUFFER_CAP - file_buf_cap) {
    memcpy(file_buf + file_buf_cap, next_file_buf, next_file_buf_cap);
    file_buf_cap += next_file_buf_cap;

    free(next_file_buf); next_file_buf = NULL;
    next_file_buf_cap = 0;
  } else {
  // copies some data to fill file_buf to the brim
    memcpy(file_buf + file_buf_cap, next_file_buf, FILE_BUFFER_CAP - file_buf_cap);
    file_buf_cap = FILE_BUFFER_CAP;

    // Move data around in secondary buffer so that it starts at zero:
    size_t residual2 = next_file_buf_cap - (FILE_BUFFER_CAP - residual);
    
    residual_data_buffer = (char*)malloc(residual2 + 1);
    memcpy(residual_data_buffer, next_file_buf + FILE_BUFFER_CAP - residual, residual2);
  
    next_file_buf = (char*)realloc(file_tmp = next_file_buf, FILE_BUFFER_CAP + 1);
    memcpy(next_file_buf, residual_data_buffer, residual2);
    free(residual_data_buffer);  
    next_file_buf_cap = residual2;
  }
  return(0);
}

inline int pbam_in::read_file_to_buffer(char * buf, const size_t len) {
  std::vector<size_t> len_chunks;
  std::vector<size_t> len_starts;
  
  if(multiFileRead && threads_to_use > 1 && FILENAME.size() > 0) {
    // Assign starts and lengths for each thread   
    size_t cur_len = 0;
    size_t chunk_div = (len+1) / threads_to_use;
    for(unsigned int i = 0; i < threads_to_use; i++) {
      // Rcpp::Rcout << "thread " << i << ", file pos " << cur_len << '\n';
      len_starts.push_back( cur_len );
      len_chunks.push_back( std::min( len - cur_len, chunk_div) );
      cur_len += len_chunks.at(i);
    }
    if(cur_len < len) {
      len_chunks.at(threads_to_use - 1) = len - len_starts.at(threads_to_use - 1);
    }
    // Spawn child ifstreams
    std::vector<std::ifstream> INchild(threads_to_use);
    const size_t cur_pos = (size_t)tellg();
    IN->seekg(len, std::ios_base::cur);

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(threads_to_use) schedule(static,1)
    #endif
    for(unsigned int k = 0; k < threads_to_use; k++) {
      INchild.at(k).open(FILENAME, std::ios::in | std::ifstream::binary);
      INchild.at(k).seekg(cur_pos + len_starts.at(k), std::ios_base::beg);
      INchild.at(k).read(buf + len_starts.at(k), len_chunks.at(k));
      INchild.at(k).close();
    }
  } else {
    IN->read(buf, len);
  }
  return(0);
}

inline size_t pbam_in::load_from_file(const size_t n_bytes) {
  /*  
    Read from file to fill n_bytes of file_buf
    First removes data upstream of file cursor
    Then fills buffer to n_bytes
  */

  char * file_tmp;
  char * residual_data_buffer;
  size_t residual = file_buf_cap-file_buf_cursor;

  size_t n_bytes_to_load =  std::min( std::max(n_bytes, residual) , FILE_BUFFER_CAP);  // Cap at file buffer
  size_t n_bytes_to_read = std::min(n_bytes_to_load - residual, IS_LENGTH - tellg());
  if(n_bytes_to_read == 0) return(0);
  // Rcpp::Rcout << "\nload_from file: n_bytes_to_read = " << n_bytes_to_read << '\n';
  // Remove residual bytes to beginning of buffer:
  if(residual > 0) {
    // have to move bytes around
    residual_data_buffer = (char*)malloc(residual + 1);
    memcpy(residual_data_buffer, file_buf + file_buf_cursor, residual);
  
    file_buf = (char*)realloc(file_tmp = file_buf, n_bytes_to_read + residual + 1);
    memcpy(file_buf, residual_data_buffer, residual);
    free(residual_data_buffer);
  } else {
    // residual == 0
    file_buf = (char*)realloc(file_tmp = file_buf, n_bytes_to_read + 1);
  }
  file_buf_cursor = 0;
  file_buf_cap = residual;
  
  // IN->read(file_buf + file_buf_cap, n_bytes_to_read);
  read_file_to_buffer(file_buf + file_buf_cap, n_bytes_to_read);

  file_buf_cap += n_bytes_to_read;
  return(n_bytes_to_read);
}

inline size_t pbam_in::fill_file_buffer() {
  return(load_from_file(FILE_BUFFER_CAP));
}

inline size_t pbam_in::read_file_chunk_to_spare_buffer(const size_t n_bytes) {
  // Read from file to fill n_bytes of next_file_buf
  if(n_bytes <= next_file_buf_cap) return(0);
  if(FILE_BUFFER_CAP <= next_file_buf_cap) return(0);
  
  size_t n_bytes_to_load =  std::min( std::max(n_bytes, next_file_buf_cap) , FILE_BUFFER_CAP);  // Cap at file buffer
  size_t n_bytes_to_read = std::min(n_bytes_to_load - next_file_buf_cap, IS_LENGTH - tellg());  

  if(n_bytes_to_read == 0) return(0);
  
  char * file_tmp;
  if(next_file_buf_cap == 0) {
    next_file_buf = (char*)realloc(file_tmp = next_file_buf, FILE_BUFFER_CAP + 1);
  }
  
  // IN->read(next_file_buf + next_file_buf_cap, n_bytes_to_read);
  read_file_to_buffer(next_file_buf + next_file_buf_cap, n_bytes_to_read);
  
  next_file_buf_cap += n_bytes_to_read;
  return(n_bytes_to_read);
}

inline int pbam_in::clean_data_buffer(const size_t n_bytes_to_decompress) {
  // Remove residual bytes to beginning of buffer:
  // Rcpp::Rcout << "clean_data_buffer\n";
  char * data_tmp;
  char * temp_buffer;
  size_t residual = data_buf_cap-data_buf_cursor;
  
  if(residual > 0) {
    // have to move bytes around
    temp_buffer = (char*)malloc(residual + 1);
    memcpy(temp_buffer, data_buf + data_buf_cursor, residual);
  
    data_buf = (char*)realloc(data_tmp = data_buf, std::max(n_bytes_to_decompress + 1, residual + 1));
    memcpy(data_buf, temp_buffer, residual);
    free(temp_buffer);
    
    data_buf_cap = residual;
  } else {
    data_buf = (char*)realloc(data_tmp = data_buf, n_bytes_to_decompress + 1);
    data_buf_cap = 0;
  }
  data_buf_cursor = 0;

  return(0);
}

inline size_t pbam_in::decompress(const size_t n_bytes_to_decompress) {
/*
  Wipes all data behind data_buf_cursor
  If file_buf data remaining is below residual (FILE_BUFFER_CAP / chunks_per_file_buf)
    then runs swap_file_buffer
  If file_buf_cursor / chunks_per_file_buf > next_file_buf_cap / chunks_per_file_buf
    then runs read_file_chunk_to_spare_buffer
  Decompresses any data in file_buf to fill up to n_bytes_to_decompress in data_buf
*/
  // Rcpp::Rcout << "decompress\n";
  clean_data_buffer(n_bytes_to_decompress);
  if(n_bytes_to_decompress < data_buf_cap) return(0);
  
  size_t decomp_cursor = data_buf_cap;  // The cursor to the data buffer to begin adding data
  size_t max_bytes_to_decompress = std::min(n_bytes_to_decompress, DATA_BUFFER_CAP) - decomp_cursor;
  // Rcpp::Rcout << "data_buf_cap " << data_buf_cap << " n_bytes_to_decompress " << n_bytes_to_decompress << "\n";

// Check if file buffer needs filling
  size_t spare_bytes_to_fill = 0;
  size_t chunk_size = (size_t)(FILE_BUFFER_CAP / chunks_per_file_buf);
  if(!eof()) {
    // If asking to decompress a large amount of data:
    if(next_file_buf_cap == 0) {
      // If secondary buffer is not yet in play, fill primary buffer first:
      if(max_bytes_to_decompress > chunk_size) {
        // Obviously asking to decompress a lot of data. Use full buffer
        load_from_file(FILE_BUFFER_CAP);
        spare_bytes_to_fill = std::min(chunk_size, IS_LENGTH - (size_t)tellg());
      } else {
        // If only asking for small amount of data, do not use full buffer
        // This is typically only called when header is read
        load_from_file(max_bytes_to_decompress);
      }
    } else if(next_file_buf_cap > 0) {
      // Checks if need to swap file buffer
      swap_file_buffer_if_needed();
      if(file_buf_cursor + chunk_size > file_buf_cap - chunk_size) {
        spare_bytes_to_fill = FILE_BUFFER_CAP;
      } else {
        spare_bytes_to_fill = file_buf_cursor + chunk_size;
      }
    }
  } else {
    // If EOF, only check when buffer needs to be swapped
    swap_file_buffer_if_needed();
    spare_bytes_to_fill = 0;
    if(file_buf_cap == file_buf_cursor) return(0);
  }

  // Set decomp_threads = threads_to_use - 1 to use asynchronous file reading
  unsigned int decomp_threads = threads_to_use;
  if(spare_bytes_to_fill > 0) {
    if(decomp_threads > 1 && !multiFileRead) {
        decomp_threads--;
    } else {
      read_file_chunk_to_spare_buffer(spare_bytes_to_fill);
      spare_bytes_to_fill = 0;
    }
  }

  size_t src_cursor = 0;
  size_t dest_cursor = 0;
  
  std::vector<size_t> src_bgzf_pos;   // source cursors
  std::vector<size_t> dest_bgzf_pos;  // dest cursors    
  std::vector<size_t> src_bgzf_cap;   // source cap
  std::vector<size_t> dest_bgzf_cap;  // dest cap    
  bool check_gzip_head = false;
  uint32_t * u32; uint16_t * u16;

  // Trial run to see if data limited at source file cap, or destination cap
  unsigned int bgzf_count = 0;
  while(src_cursor < chunk_size) {
    if(file_buf_cursor + src_cursor + 28 > file_buf_cap) break;

    // Check bamGzipHead every 1000 blocks
    if(bgzf_count % 1000 == 0) check_gzip_head = true; // Check every 10000 blocks
    bgzf_count++;
    
    if(check_gzip_head) {
      // Rcpp::Rcout << "file_buf_cursor " << file_buf_cursor << " src_cursor " << src_cursor << '\n';
      if(strncmp(bamGzipHead, file_buf + file_buf_cursor + src_cursor, bamGzipHeadLength) != 0) {
        break;  // Gzip head corrupt
      } else {
        check_gzip_head = false;
      }
    }

    u16 = (uint16_t*)(file_buf + file_buf_cursor + src_cursor + 16);
    if(file_buf_cursor + src_cursor + (*u16+1) > file_buf_cap) break;

    if(src_cursor + (*u16+1) > chunk_size) break; // will not read above chunk_size
    
    // Abort if filling dest with extra data will lead to overflow
    u32 = (uint32_t*)(file_buf + file_buf_cursor + src_cursor + (*u16+1) - 4);
    if(dest_cursor + *u32 > max_bytes_to_decompress) break;   
    
    // Rcpp::Rcout << "src_cursor = " << src_cursor << " dest_cursor = " << dest_cursor << '\n';
    src_cursor += *u16 + 1;
    dest_cursor += *u32;

  }
  // Rcpp::Rcout << "src_cursor = " << src_cursor << " dest_cursor = " << dest_cursor << "done\n";

  if(check_gzip_head) {
    Rcpp::Rcout << "BGZF blocks corrupt\n";
    return(0);
  }

// Usage of source bytes to divide workload:
  size_t src_max = src_cursor; size_t dest_max = dest_cursor;
  size_t src_divider = 1 + (src_max / decomp_threads);
  size_t next_divider = src_divider;
  src_cursor = 0;
  dest_cursor = 0;
  
  // Profile the file and increment until the destination buffer is reached

  // Rcpp::Rcout << "file_buf_cursor = " << file_buf_cursor <<
    // " file_buf_cap = " << file_buf_cap <<
    // " next_file_buf_cap = " << next_file_buf_cap <<
    // " data_buf_cap = " << data_buf_cap << 
    // " data_buf_cursor = " << data_buf_cursor << '\n';
  unsigned int threads_accounted_for = 0;
  src_bgzf_pos.push_back(file_buf_cursor + src_cursor);  
  dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
  while(src_cursor < src_max) {
    u16 = (uint16_t*)(file_buf + file_buf_cursor + src_cursor + 16);
    u32 = (uint32_t*)(file_buf + file_buf_cursor + src_cursor + (*u16+1) - 4);
    src_cursor += *u16 + 1;
    dest_cursor += *u32;
    

    if(src_cursor > next_divider && threads_accounted_for < decomp_threads - 1) {
      src_bgzf_cap.push_back(file_buf_cursor + src_cursor);  
      dest_bgzf_cap.push_back(decomp_cursor + dest_cursor);
      src_bgzf_pos.push_back(file_buf_cursor + src_cursor);  
      dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
      next_divider = std::min(next_divider + src_divider, src_max);
      threads_accounted_for++;
      // Rcpp::Rcout << "src_cursor = " << src_cursor << " dest_cursor = " << dest_cursor << '\n';
    }
  }
  
  while(threads_accounted_for < decomp_threads - 1) {
    src_bgzf_cap.push_back(file_buf_cursor + src_cursor);  
    dest_bgzf_cap.push_back(decomp_cursor + dest_cursor);
    src_bgzf_pos.push_back(file_buf_cursor + src_cursor);  
    dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
    // next_divider = std::min(next_divider + src_divider, src_max);
    threads_accounted_for++;
    // Rcpp::Rcout << "src_cursor = " << src_cursor << " dest_cursor = " << dest_cursor << '\n';
  }  
  // Rcpp::Rcout << "src_cursor = " << src_cursor << " dest_cursor = " << dest_cursor << "done\n";

  // NB: last src and dest bgzf pos are the end: So always one more bgzf pos than done
  src_bgzf_cap.push_back(file_buf_cursor + src_max);   
  dest_bgzf_cap.push_back(decomp_cursor + dest_max);
   
  // Now comes the multi-threaded decompression:
  bool error_occurred = false;
  
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads_to_use) schedule(static,1)
  #endif
  for(unsigned int k = 0; k < threads_to_use; k++) {
    if(k == decomp_threads) {
      // as decomp_threads == threads_to_use - 1 iff spare_bytes_to_fill > 0
      
      // Last thread is used to prime spare file buffer
      read_file_chunk_to_spare_buffer(spare_bytes_to_fill);
      spare_bytes_to_fill = 0;
    } else {
      // unsigned int i = thrd_partitions.at(k);

      size_t thread_src_cursor = src_bgzf_pos.at(k);
      size_t thread_dest_cursor = dest_bgzf_pos.at(k);
      uint32_t crc = 0;
      uint32_t * crc_check;
      uint16_t * src_size;
      uint32_t * dest_size;

      z_stream * zs;
      while(thread_src_cursor < src_bgzf_cap.at(k) && !error_occurred) {
        src_size = (uint16_t *)(file_buf + thread_src_cursor + 16);
        crc_check = (uint32_t *)(file_buf + thread_src_cursor + *src_size+1 - 8);
        dest_size = (uint32_t *)(file_buf + thread_src_cursor + *src_size+1 - 4);

        if(*dest_size > 0) {
          // Rcpp::Rcout << thread_src_cursor << "\t" << thread_dest_cursor << '\t' << *src_size+1 << "\t" << *dest_size << '\n';
          zs = new z_stream;
          zs->zalloc = NULL; zs->zfree = NULL; zs->msg = NULL;
          zs->next_in = (Bytef*)(file_buf + thread_src_cursor + 18);
          zs->avail_in = *src_size + 1 - 18;
          zs->next_out = (Bytef*)(data_buf + thread_dest_cursor);
          zs->avail_out = *dest_size;

          int ret = inflateInit2(zs, -15);
          if(ret != Z_OK) {
            Rcpp::Rcout << "Exception during BAM decompression - inflateInit2() fail: (" << ret << ") \n";
              
  #ifdef _OPENMP
  #pragma omp critical
  #endif
            error_occurred = true;
          }
          if(!error_occurred) {
            ret = inflate(zs, Z_FINISH);
            if(ret != Z_OK && ret != Z_STREAM_END) {
              Rcpp::Rcout << "Exception during BAM decompression - inflate() fail: (" << ret << ") \n";
  #ifdef _OPENMP
  #pragma omp critical
  #endif
              error_occurred = true;
            }
          }
          if(!error_occurred) {
            ret = inflateEnd(zs);

            crc = crc32(crc32(0L, NULL, 0L), 
            (Bytef*)(data_buf + thread_dest_cursor), *dest_size);
            if(*crc_check != crc) {
              Rcpp::Rcout << "CRC fail during BAM decompression\n";
  #ifdef _OPENMP
  #pragma omp critical
  #endif
              error_occurred = true;
            }

          }
          
          delete zs;
        }
        thread_src_cursor += *src_size + 1;
        thread_dest_cursor += *dest_size;
// Rcpp::Rcout << "thread_src_cursor = " << thread_src_cursor << " thread_dest_cursor = " << thread_dest_cursor << '\n';

      }
    }
    
  }
  
  if(error_occurred) {
    Rcpp::Rcout << "error_occurred\n";
    return(0);
  }
  
  file_buf_cursor = src_bgzf_cap.at(src_bgzf_cap.size() - 1);
  data_buf_cap = dest_bgzf_cap.at(dest_bgzf_cap.size() - 1);
  // Rcpp::Rcout << "file_buf_cursor = " << file_buf_cursor << " data_buf_cap = " << data_buf_cap << '\n';
  // Rcpp::Rcout << "decompress() done\ttellg() = " << tellg() << '\n';
  return(dest_cursor);
}

inline unsigned int pbam_in::read(char * dest, const unsigned int len) {
  // Reads bytes from the current data_buf, increments data_buff_cursor
  if(data_buf_cap - data_buf_cursor < len) {
    decompress(len + 65536);
  }
  
  unsigned int n_bytes_to_read = std::min((size_t)len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  // Rcpp::Rcout << n_bytes_to_read << '\n';
  memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}

inline unsigned int pbam_in::ignore(const unsigned int len) {
  // Ignores bytes from the current data_buf, increments data_buff_cursor
  if(data_buf_cap - data_buf_cursor < len) {
    decompress(len + 65536);
  }
  
  unsigned int n_bytes_to_read = std::min((size_t)len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  // memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}

inline unsigned int pbam_in::peek(char * dest, const unsigned int len) {
  // Peeks at bytes from the current data_buf, does not increment data_buff_cursor
  if(data_buf_cap - data_buf_cursor < len) {
    decompress(len + 65536);
  }
  
  unsigned int n_bytes_to_read = std::min((size_t)len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  // data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
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
  // Rcpp::Rcout << " n_ref " << n_ref << '\n';
  char chrom_buffer[1000]; std::string chrName;
  for(unsigned int i = 0; i < n_ref; i++) {
    read(u32c, 4); u32 = (uint32_t *)u32c;
    read(chrom_buffer, *u32);
    chrName = std::string(chrom_buffer, *u32-1);
    chr_names.push_back(chrName);
    // Rcpp::Rcout << " chr_names " << chrName << '\n';
    read(u32c, 4); u32 = (uint32_t *)u32c;
    chr_lens.push_back(*u32);
  }
  
  free(u32c);
  return(0);
}

inline int pbam_in::obtainChrs(std::vector<std::string> & s_chr_names, std::vector<uint32_t> & u32_chr_lens) {
  if(!magic_header) {
    Rcpp::Rcout << "Header is not yet read\n";
    return(-1);
  }
  if(n_ref == 0) {
    Rcpp::Rcout << "No chromosome names stored. Is pbam_in::readHeader() been run yet?\n";
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

inline int pbam_in::fillReads() {
  // Returns -1 if error, and 1 if EOF. Otherwise, returns 0
  // Rcpp::Rcout << "fillReads()\n";
  // If header not read:
  if(!magic_header) {
    Rcpp::Rcout << "Header is not yet read\n";
    return(-1);
  }
  if(n_ref == 0) {
    Rcpp::Rcout << "No chromosome names stored. Is pbam_in::readHeader() been run yet?\n";
    return(-1);
  }
  
  // Check if previous reads are all read:
  if(read_cursors.size() > 0) {
    for(unsigned int i = 0; i < read_cursors.size(); i++) {
      if(read_cursors.at(i) < read_ptr_ends.at(i)) {
        Rcpp::Rcout << "Thread " << i << " has reads remaining. Please debug your code ";
        Rcpp::Rcout << "and make sure all threads clear their reads before filling any more reads\n";
        return(-1);
      }
    }
  }
  
  // Clear read pointers:
  read_cursors.resize(0);
  read_ptr_ends.resize(0);
  
  decompress(DATA_BUFFER_CAP);
  // Rcpp::Rcout << "Decompress done\n";
  
  // Check decompressed data contains at least 1 full read
  uint32_t *u32p;
  if(data_buf_cap - data_buf_cursor < 4) {
    return(1);
  } else {
    u32p = (uint32_t *)(data_buf + data_buf_cursor);
    if(*u32p + 4 > data_buf_cap - data_buf_cursor) {
      return(1);
    }
  }
  
  // Roughly divide the buffer into N regions:
  size_t data_divider = 1 + ((data_buf_cap - data_buf_cursor) / threads_to_use);
  size_t next_divider = std::min(data_buf_cursor + data_divider, data_buf_cap);
  // Iterates through data buffer aand assigns pointers to beginning of reads
  // Rcpp::Rcout << "data_buf_cursor = " << data_buf_cursor << " data_buf_cap = " << data_buf_cap << '\n';
  // Rcpp::Rcout << "tellg() = " << tellg() << " file_buf_cap = " << file_buf_cap << 
    // " file_buf_cursor = " << file_buf_cursor << " next_file_buf_cap = " << next_file_buf_cap << '\n';
  // bool has_reads_left_in_buffer = true;
  read_cursors.push_back(data_buf_cursor);
  unsigned int threads_accounted_for = 0;
  while(1) {
    // Checks remaining data contains at least 1 full read; breaks otherwise
    if(data_buf_cap - data_buf_cursor >= 4) {
      u32p = (uint32_t *)(data_buf + data_buf_cursor);
      if(*u32p + 4 <= data_buf_cap - data_buf_cursor) {
        data_buf_cursor += *u32p + 4;
      } else {
        break;
      }
    } else {
      break;
    }
    if(data_buf_cursor >= next_divider && threads_accounted_for < threads_to_use - 1) {
      read_ptr_ends.push_back(data_buf_cursor);
      read_cursors.push_back(data_buf_cursor);
      next_divider = std::min(data_buf_cap, next_divider + data_divider);
      threads_accounted_for++;
      // Rcpp::Rcout << data_buf_cursor << '\t';
    }
  }

  while(threads_accounted_for < threads_to_use - 1) {
    read_ptr_ends.push_back(data_buf_cursor);
    read_cursors.push_back(data_buf_cursor);
    // next_divider = std::max(data_buf_cap, next_divider + data_divider);
    threads_accounted_for++;
    // Rcpp::Rcout << data_buf_cursor << '\t';
  }

  read_ptr_ends.push_back(data_buf_cursor);
  // Rcpp::Rcout << data_buf_cursor << '\n';
  return(0);
}

inline pbam1_t pbam_in::supplyRead(const unsigned int thread_number) {
  pbam1_t read;
  if(thread_number > read_cursors.size()) {
    Rcpp::Rcout << "Invalid thread number parsed to supplyRead()\n";
    return(read);
  }
  if(read_cursors.at(thread_number) >= read_ptr_ends.at(thread_number)) {
    return(read);
  }
  read = pbam1_t(data_buf + read_cursors.at(thread_number), false);
  if(read.validate()) {
    read_cursors.at(thread_number) += read.block_size() + 4;
  } else {
    // Check this is actually end of thread buffer; throw error here otherwise
    if(read_cursors.at(thread_number) < read_ptr_ends.at(thread_number)) {
      Rcpp::Rcout << "Invalid read found before end of thread buffer " 
        << thread_number << ". read_cursor = " << read_cursors.at(thread_number)
        << ", read_ptr_ends = " << read_ptr_ends.at(thread_number) << '\n';
    }
  }    
  return(read);
}

#endif