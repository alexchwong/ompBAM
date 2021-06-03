#ifndef _pbam_in
#define _pbam_in

class pbam_in {
  private:
// pbam_in settings
    size_t FILE_BUFFER_CAP = 200000000;     // 1e8
    size_t DATA_BUFFER_CAP = 1000000000;    // 1e9
    unsigned int FILE_BUFFER_SEGMENTS = 10;  // Divide file buffer into n segments
    unsigned int threads_to_use = 1;
  
// File particulars
    std::istream * IN;    
    size_t IS_LENGTH;     // size of BAM file (Input Stream Length)
  
    // Header storage:
    char * magic_header;  // Always of size = 8
    uint32_t l_text;      // sizeof *headertext
    char * headertext;
    uint32_t n_ref;
    std::vector<std::string> chr_names;
    std::vector<uint32_t> chr_lens;
  
    // Data buffers
    char * file_buf; size_t file_buf_cap; size_t file_buf_cursor;
    char * next_file_buf; size_t next_file_buf_cap;
    
    char * data_buf; size_t data_buf_cap; size_t data_buf_cursor;
    
    std::vector<char *> read_ptrs;              // Raw char pointers to data_buf
    std::vector<uint32_t> read_ptr_partitions;  // std::maximum read # for threads 0 to i-1
    std::vector<uint32_t> read_cursors;         // Read # to be returned @ next call to supplyRead

    void initialize_buffers();
    void clear_buffers();

    int swap_file_buffer_if_needed();
    size_t load_from_file(size_t n_bytes);
    size_t fill_file_buffer();

    unsigned int calculate_chunks_to_load_to_secondary_buffer();
    size_t read_file_chunk_to_spare_buffer(size_t n_bytes);

    size_t decompress(size_t n_bytes_to_decompress);
    
    
    int clean_data_buffer(size_t n_bytes_to_decompress);

    // Only for reading header
    unsigned int read(char * dest, unsigned int len);  // returns the number of bytes actually read
    unsigned int ignore(unsigned int len);
    unsigned int peek(char * dest, unsigned int len) ;

    // Convenience functions
    size_t tellg() {return((size_t)IN->tellg());};
    bool eof() {return(IS_LENGTH == tellg());};
    bool fail() {return(IN->fail());};
    
    // Disable copy construction / assignment:
    pbam_in(const pbam_in &t);
    pbam_in & operator = (const pbam_in &t);
  public:
    pbam_in();
    pbam_in(size_t file_buffer_cap, size_t data_buffer_cap, unsigned int file_buffer_segments);
    ~pbam_in();

    int SetInputHandle(std::istream *in_stream, unsigned int n_threads); // opens file stream, checks EOF and file size
    int readHeader();
    int obtainChrs(std::vector<std::string> & s_chr_names, std::vector<uint32_t> & u32_chr_lens);

    size_t GetLength() { return(IS_LENGTH); };

    int fillReads();  // Returns 1 if no more reads available
    
    // Read functions:
    pbam1_t supplyRead(unsigned int thread_number = 0);
};

// Functions

void pbam_in::initialize_buffers() {
  file_buf = NULL;
  data_buf = NULL;
  next_file_buf = NULL;
  file_buf_cap = 0; file_buf_cursor = 0;
  data_buf_cap = 0; data_buf_cursor = 0;
  next_file_buf_cap = 0;
  magic_header = NULL; l_text = 0; headertext = NULL; n_ref = 0;
  chr_names.resize(0); chr_lens.resize(0);

  read_ptrs.resize(0);
  read_ptr_partitions.resize(0);
  read_cursors.resize(0);
  
  IN = NULL;
}

void pbam_in::clear_buffers() {
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
  
  next_file_buf = NULL;
  file_buf_cap = 0; file_buf_cursor = 0;
  data_buf_cap = 0; data_buf_cursor = 0;
  next_file_buf_cap = 0;
  magic_header = NULL; l_text = 0; headertext = NULL; n_ref = 0;
  chr_names.resize(0); chr_lens.resize(0);

  read_ptrs.resize(0);
  read_ptr_partitions.resize(0);
  read_cursors.resize(0);
  
  IN = NULL;
}

pbam_in::pbam_in() {
  // Initialize buffers
  initialize_buffers();
  // Settings reset
  // FILE_BUFFER_CAP = 1000000000;
  // DATA_BUFFER_CAP = 1000000000;
  // FILE_BUFFER_SEGMENTS = 5;
  // threads_to_use = 1;
  
  IN = NULL;
}

pbam_in::pbam_in(size_t file_buffer_cap, size_t data_buffer_cap, unsigned int file_buffer_segments) {
  // Initialize buffers
  initialize_buffers();
  
  FILE_BUFFER_CAP = file_buffer_cap;
  DATA_BUFFER_CAP = data_buffer_cap;
  FILE_BUFFER_SEGMENTS = file_buffer_segments;
  threads_to_use = 1;
  
  IN = NULL;  
}

pbam_in::~pbam_in() {
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
}

int pbam_in::SetInputHandle(std::istream *in_stream, unsigned int n_threads) {
  #ifdef _OPENMP
    if(n_threads > (unsigned int)omp_get_max_threads()) {
      threads_to_use = (unsigned int)omp_get_max_threads();
    } else {
      threads_to_use = n_threads;
    }
  #else
    threads_to_use = 1;
  #endif
  
  IN = in_stream;
  
  // Assign IS_LENGTH
  IN->seekg(0, std::ios_base::end);
  IS_LENGTH = IN->tellg();
  
  // Check valid EOF:
  IN->seekg(IS_LENGTH-bamEOFlength, std::ios_base::beg);
  
  char eof_check[bamEOFlength + 1];
  IN->read(eof_check, bamEOFlength);
  if(strncmp(bamEOF, eof_check, bamEOFlength) != 0) {
    Rcout << "Error opening BAM - EOF bit corrupt. Perhaps this file is truncated?\n";
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
}

int pbam_in::swap_file_buffer_if_needed() {
  // Transfers residual data from current file buffer to the next:
  Rcout << "swap_file_buffer_if_needed()\n";
  if(next_file_buf_cap == 0) return(1);
  size_t chunk_size = (size_t)(FILE_BUFFER_CAP / FILE_BUFFER_SEGMENTS);
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
  
  if(next_file_buf_cap <= FILE_BUFFER_CAP - file_buf_cap) {
    memcpy(file_buf + file_buf_cap, next_file_buf, next_file_buf_cap);
    file_buf_cap += next_file_buf_cap;

    free(next_file_buf); next_file_buf = NULL;
    next_file_buf_cap = 0;
  } else {
    memcpy(file_buf + file_buf_cap, next_file_buf, FILE_BUFFER_CAP - file_buf_cap);
    file_buf_cap = FILE_BUFFER_CAP;

    // Move data around in secondary buffer:
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


size_t pbam_in::load_from_file(size_t n_bytes) {
  /*  
    Read from file to fill n_bytes of file_buf
    First removes data upstream of file cursor
    Then fills buffer to n_bytes
  */
  // Rcout << "load from file()\n";
  char * file_tmp;
  char * residual_data_buffer;
  size_t residual = file_buf_cap-file_buf_cursor;
  // Rcout << "Residual = " << residual << '\n';
  size_t n_bytes_to_load =  std::min( std::max(n_bytes, residual) , FILE_BUFFER_CAP);  // Cap at file buffer
  size_t n_bytes_to_read = std::min(n_bytes_to_load - residual, IS_LENGTH - IN->tellg());
  if(n_bytes_to_read == 0) return(0);
  
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
  
  IN->read(file_buf + file_buf_cap, n_bytes_to_read);
  file_buf_cap += n_bytes_to_read;
  // Rcout << "load from file() finished - n_bytes_to_read = " << n_bytes_to_read << "\n";
  return(n_bytes_to_read);
}

size_t pbam_in::fill_file_buffer() {
  return(load_from_file(FILE_BUFFER_CAP));
}

unsigned int pbam_in::calculate_chunks_to_load_to_secondary_buffer() {
  size_t chunk_size = (size_t)(FILE_BUFFER_CAP / FILE_BUFFER_SEGMENTS);
  unsigned int n_chunks_to_fetch = (file_buf_cursor / chunk_size) + 1;
  unsigned int n_chunks_avail = (next_file_buf_cap / chunk_size);
  if(n_chunks_to_fetch <= n_chunks_avail) return(0);
  return(n_chunks_to_fetch - n_chunks_avail);
}

size_t pbam_in::read_file_chunk_to_spare_buffer(size_t n_bytes) {
  // Read from file to fill n_bytes of next_file_buf
  if(n_bytes <= next_file_buf_cap) return(0);
  if(FILE_BUFFER_CAP <= next_file_buf_cap) return(0);
  
  size_t n_bytes_to_load =  std::min( std::max(n_bytes, next_file_buf_cap) , FILE_BUFFER_CAP);  // Cap at file buffer
  size_t n_bytes_to_read = std::min(n_bytes_to_load - next_file_buf_cap, IS_LENGTH - IN->tellg());  

  if(n_bytes_to_read == 0) return(0);
  
  char * file_tmp;
  if(next_file_buf_cap == 0) {
    next_file_buf = (char*)realloc(file_tmp = next_file_buf, FILE_BUFFER_CAP + 1);
  }
  
  IN->read(next_file_buf + next_file_buf_cap, n_bytes_to_read);
  next_file_buf_cap += n_bytes_to_read;
  return(n_bytes_to_read);
}

int pbam_in::clean_data_buffer(size_t n_bytes_to_decompress) {
  // Remove residual bytes to beginning of buffer:
  // Rcout << "clean_data_buffer\n";
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

size_t pbam_in::decompress(size_t n_bytes_to_decompress) {
/*
  Wipes all data behind data_buf_cursor
  If file_buf data remaining is below residual (FILE_BUFFER_CAP / FILE_BUFFER_SEGMENTS)
    then runs swap_file_buffer
  If file_buf_cursor / FILE_BUFFER_SEGMENTS > next_file_buf_cap / FILE_BUFFER_SEGMENTS
    then runs read_file_chunk_to_spare_buffer
  Decompresses any data in file_buf to fill up to n_bytes_to_decompress in data_buf
*/
  // Rcout << "decompress\n";
  clean_data_buffer(n_bytes_to_decompress);
  if(n_bytes_to_decompress < data_buf_cap) return(0);
  
  size_t decomp_cursor = data_buf_cap;  // The cursor to the data buffer to begin adding data
  size_t max_bytes_to_decompress = std::min(n_bytes_to_decompress, DATA_BUFFER_CAP) - decomp_cursor;
  // Rcout << "data_buf_cap " << data_buf_cap << " n_bytes_to_decompress " << n_bytes_to_decompress << "\n";

// Check if file buffer needs filling
  unsigned int n_chunks_to_load = 0;
  size_t spare_bytes_to_fill = 0;
  size_t chunk_size = (size_t)(FILE_BUFFER_CAP / FILE_BUFFER_SEGMENTS);
  if(!eof()) {
    // If asking to decompress a large amount of data:
    if(next_file_buf_cap == 0) {
      // If secondary buffer is not yet in play, fill primary buffer first:
      load_from_file(max_bytes_to_decompress);
      if(max_bytes_to_decompress > FILE_BUFFER_CAP) {
        // Ask to load secondary buffer from file
        // Note that spare buffer is never filled unless bytes_to_decompress > FILE_BUFFER_CAP
        n_chunks_to_load = calculate_chunks_to_load_to_secondary_buffer();
      }
    } else if(next_file_buf_cap > 0) {
      // Checks if need to swap file buffer
      swap_file_buffer_if_needed();
      n_chunks_to_load = calculate_chunks_to_load_to_secondary_buffer();
    }
    if(n_chunks_to_load > 0) {
      unsigned int n_chunks_avail = (next_file_buf_cap / chunk_size);
      spare_bytes_to_fill = chunk_size * (n_chunks_avail + n_chunks_to_load);
    }
  }

  size_t src_cursor = 0;
  size_t dest_cursor = 0;
  
  std::vector<size_t> src_bgzf_pos;
  std::vector<size_t> dest_bgzf_pos;
  bool check_gzip_head = false;
  uint32_t * u32; uint16_t * u16;

  // Profile the file and increment until the destination buffer is reached
  unsigned int bgzf_count = 0;
  Rcout << "file_buf_cursor = " << file_buf_cursor <<
    " file_buf_cap = " << file_buf_cap <<
    " next_file_buf_cap = " << next_file_buf_cap <<
    " data_buf_cap = " << data_buf_cap << 
    " data_buf_cursor = " << data_buf_cursor << '\n';
  while(1) {
    // Abort if cannot read smallest possible bgzf block:
    if(file_buf_cursor + src_cursor + 28 > file_buf_cap) break;

    // Check bamGzipHead every 1000 blocks
    if(bgzf_count % 1000 == 0) check_gzip_head = true; // Check every 10000 blocks
    if(check_gzip_head) {
      // Rcout << "file_buf_cursor " << file_buf_cursor << " src_cursor " << src_cursor << '\n';
      if(strncmp(bamGzipHead, file_buf + file_buf_cursor + src_cursor, bamGzipHeadLength) != 0) {
        break;  // Gzip head corrupt
      } else {
        check_gzip_head = false;
      }
    }
    
    // Abort if cannot read entire bgzf block:
    u16 = (uint16_t*)(file_buf + file_buf_cursor + src_cursor + 16);
    if(file_buf_cursor + src_cursor + (*u16+1) > file_buf_cap) break;

    // Abort if filling dest with extra data will lead to overflow
    u32 = (uint32_t*)(file_buf + file_buf_cursor + src_cursor + (*u16+1) - 4);
    if(dest_cursor + *u32 > max_bytes_to_decompress) break;

    // If checks passed, designate these blocks to be decompressed
    src_bgzf_pos.push_back(file_buf_cursor + src_cursor);  
    src_cursor += *u16 + 1;
    
    dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
    dest_cursor += *u32;

    bgzf_count++;
  }

  if(check_gzip_head) {
    Rcout << "BGZF blocks corrupt\n";
    return(0);
  }
  
  // NB: last src and dest bgzf pos are the end: So always one more bgzf pos than done
  src_bgzf_pos.push_back(file_buf_cursor + src_cursor);   
  dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
  
  unsigned int decomp_threads = threads_to_use;
  if(spare_bytes_to_fill > 0) {
    if(decomp_threads > 1) {
        decomp_threads--;
    } else {
      read_file_chunk_to_spare_buffer(spare_bytes_to_fill);
      spare_bytes_to_fill = 0;
    }
  }
  
  // Assign jobs between n threads
  std::vector<unsigned int> thrd_partitions;
  unsigned int job_inc = 0;
  unsigned int job_size = (src_bgzf_pos.size()) / decomp_threads;
  for(unsigned int k = 0; k < decomp_threads; k++) {
    thrd_partitions.push_back(job_inc);
    job_inc+=job_size;
  }
  thrd_partitions.push_back(src_bgzf_pos.size()-1); 
  
  // Now comes the multi-threaded decompression:
  bool error_occurred = false;
  
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(unsigned int k = 0; k < threads_to_use; k++) {
    if(k == decomp_threads) {
      // as decomp_threads == threads_to_use - 1 iff spare_bytes_to_fill > 0
      
      // Last thread is used to prime spare file buffer
      read_file_chunk_to_spare_buffer(spare_bytes_to_fill);
      spare_bytes_to_fill = 0;
    } else {
      unsigned int i = thrd_partitions.at(k);

      size_t thread_src_cursor = 0;
      size_t thread_dest_cursor = 0;

      uint32_t crc = 0;
      uint32_t * crc_check;
      uint32_t src_size = 0;
      uint32_t dest_size = 0;

      z_stream * zs;
      while(i < thrd_partitions.at(k+1) && !error_occurred) {
        thread_src_cursor = src_bgzf_pos.at(i);
        thread_dest_cursor = dest_bgzf_pos.at(i);
        crc_check = (uint32_t *)(file_buf + src_bgzf_pos.at(i+1) - 8);
        src_size = src_bgzf_pos.at(i+1) - thread_src_cursor;
        dest_size = dest_bgzf_pos.at(i+1) - thread_dest_cursor;

        // Rcout << thread_src_cursor << "\t" << thread_dest_cursor << '\t' << src_size << "\t" << dest_size << '\n';
        zs = new z_stream;
        zs->zalloc = NULL; zs->zfree = NULL; zs->msg = NULL;
        zs->next_in = (Bytef*)(file_buf + thread_src_cursor + 18);
        zs->avail_in = src_size - 18;
        zs->next_out = (Bytef*)(data_buf + thread_dest_cursor);
        zs->avail_out = dest_size;

        int ret = inflateInit2(zs, -15);
        if(ret != Z_OK) {
          Rcout << "Exception during BAM decompression - inflateInit2() fail: (" << ret << ") "
            "BGZF block number" << i << '\n';
            
#ifdef _OPENMP
#pragma omp critical
#endif
          error_occurred = true;
        }
        if(!error_occurred) {
          ret = inflate(zs, Z_FINISH);
          if(ret != Z_OK && ret != Z_STREAM_END) {
            Rcout << "Exception during BAM decompression - inflate() fail: (" << ret << ") "
              "BGZF block number" << i << '\n';
#ifdef _OPENMP
#pragma omp critical
#endif
            error_occurred = true;
          }
        }
        if(!error_occurred) {
          ret = inflateEnd(zs);

          crc = crc32(crc32(0L, NULL, 0L), 
          (Bytef*)(data_buf + thread_dest_cursor), dest_size);
          if(*crc_check != crc) {
            Rcout << "CRC fail during BAM decompression" << ", BGZF block number" << i << '\n';
#ifdef _OPENMP
#pragma omp critical
#endif
            error_occurred = true;
          }

        }
        
        delete zs;
        i++;
      }
    }
    
  }
  
  if(error_occurred) {
    Rcout << "error_occurred\n";
    return(0);
  }
  
  file_buf_cursor = src_bgzf_pos.at(src_bgzf_pos.size() - 1);
  data_buf_cap = dest_bgzf_pos.at(dest_bgzf_pos.size() - 1);
  Rcout << "file_buf_cursor = " << file_buf_cursor << '\n';
  Rcout << "decompress() done\ttellg() = " << tellg() << '\n';
  return(dest_cursor);
}

unsigned int pbam_in::read(char * dest, unsigned int len) {
  // Reads bytes from the current data_buf, increments data_buff_cursor
  if(data_buf_cap - data_buf_cursor < len) {
    decompress(len + 65536);
  }
  
  unsigned int n_bytes_to_read = std::min((size_t)len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  // Rcout << n_bytes_to_read << '\n';
  memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}

unsigned int pbam_in::ignore(unsigned int len) {
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

unsigned int pbam_in::peek(char * dest, unsigned int len) {
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

int pbam_in::readHeader() {
  if(magic_header) {
    Rcout << "Header is already read\n";
    return(-1);
  }
  
  magic_header = (char*)malloc(8+1);
  read(magic_header, 8);
  if(strncmp(magic_header, magicstring, 4) != 0) {
    Rcout << "Invalid BAM magic string\n";
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
  // Rcout << " n_ref " << n_ref << '\n';
  char chrom_buffer[1000]; std::string chrName;
  for(unsigned int i = 0; i < n_ref; i++) {
    read(u32c, 4); u32 = (uint32_t *)u32c;
    read(chrom_buffer, *u32);
    chrName = std::string(chrom_buffer, *u32-1);
    chr_names.push_back(chrName);
    // Rcout << " chr_names " << chrName << '\n';
    read(u32c, 4); u32 = (uint32_t *)u32c;
    chr_lens.push_back(*u32);
  }
  
  free(u32c);
  return(0);
}

int pbam_in::obtainChrs(std::vector<std::string> & s_chr_names, std::vector<uint32_t> & u32_chr_lens) {
  if(!magic_header) {
    Rcout << "Header is not yet read\n";
    return(-1);
  }
  if(n_ref == 0) {
    Rcout << "No chromosome names stored. Is pbam_in::readHeader() been run yet?\n";
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

int pbam_in::fillReads() {
  // Returns -1 if error, and 1 if EOF. Otherwise, returns 0

  // If header not read:
  if(!magic_header) {
    Rcout << "Header is not yet read\n";
    return(-1);
  }
  if(n_ref == 0) {
    Rcout << "No chromosome names stored. Is pbam_in::readHeader() been run yet?\n";
    return(-1);
  }
  
  // Check if previous reads are all read:
  if(read_cursors.size() > 0) {
    for(unsigned int i = 0; i < read_cursors.size(); i++) {
      if(read_cursors.at(i) < read_ptr_partitions.at(i)) {
        Rcout << "Thread " << i << " has reads remaining. Please debug your code ";
        Rcout << "and make sure all threads clear their reads before filling any more reads\n";
        return(-1);
      }
    }
  }
  
  // Clear read pointers:
  read_ptrs.resize(0);
  read_ptr_partitions.resize(0);
  read_cursors.resize(0);

  decompress(DATA_BUFFER_CAP);
  
  // Iterates through data buffer aand assigns pointers to beginning of reads
  uint32_t *u32p;
  bool has_reads_left_in_buffer = true;
  while(has_reads_left_in_buffer) {
    u32p = (uint32_t *)(data_buf + data_buf_cursor);
    if(*u32p + 4 <= data_buf_cap - data_buf_cursor) {
      read_ptrs.push_back(data_buf + data_buf_cursor);
      data_buf_cursor += *u32p + 4;
    } else {
      has_reads_left_in_buffer = false;
    }
  }

  if(read_ptrs.size() == 0) {
    Rcout << "End of buffer reached\ttellg() = " << tellg() << '\n';
    return(1);
  }
  
  // Partition reads by number of threads
  unsigned int reads_per_thread = 1 + (read_ptrs.size() / threads_to_use);
  unsigned int cursor = 0;
  while(read_cursors.size() < threads_to_use) {
    read_cursors.push_back(cursor);
    cursor += reads_per_thread;
    if(cursor > read_ptrs.size()) cursor = read_ptrs.size();
    read_ptr_partitions.push_back(cursor);
  }
  return(0);
}

pbam1_t pbam_in::supplyRead(unsigned int thread_number) {
  pbam1_t read;
  if(thread_number > read_cursors.size()) {
    Rcout << "Invalid thread number parsed to supplyRead()\n";
    return(read);
  }
  if(read_cursors.at(thread_number) == read_ptr_partitions.at(thread_number)) {
    return(read);
  }
  read = pbam1_t(read_ptrs.at(read_cursors.at(thread_number)), false);
  read_cursors.at(thread_number)++;
  return(read);
}


#endif