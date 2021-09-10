#ifndef _pbam_in_decompress
#define _pbam_in_decompress

// All functions in this file are private functions.

/*
  Asks to decompress upto `n_bytes_to_decompress` worth of data to be contained
    in data_buf
  
  First wipes all data upstream of data_buf_cursor
  If file_buf data remaining is below residual (FILE_BUFFER_CAP / chunks_per_file_buf)
    then runs swap_file_buffer
  If file_buf_cursor / chunks_per_file_buf > next_file_buf_cap / chunks_per_file_buf
    then runs read_file_chunk_to_spare_buffer
  Decompresses any data in file_buf to fill up to n_bytes_to_decompress in data_buf
*/
inline size_t pbam_in::decompress(const size_t n_bytes_to_decompress) {
  clean_data_buffer(n_bytes_to_decompress);
  if(n_bytes_to_decompress < data_buf_cap) return(0);
  
  // The cursor to the data buffer to begin adding data
  size_t decomp_cursor = data_buf_cap;  
  size_t max_bytes_to_decompress = 
    std::min(n_bytes_to_decompress, DATA_BUFFER_CAP) - decomp_cursor;

  size_t chunk_size = (size_t)(FILE_BUFFER_CAP / chunks_per_file_buf);

  // Store an order to fill this maximum bytes to next_file_buf
  size_t spare_bytes_to_fill = 0;   
  if(!eof()) {
    if(next_file_buf_cap == 0) {
      // If secondary buffer is not yet in play, fill primary buffer first:
      if(max_bytes_to_decompress > chunk_size) {
        // If asking for more than chunk_size, fill the whole buffer
        fill_file_buffer();

        // Ask to fill chunk_size amount to next_file_buf
        spare_bytes_to_fill = std::min(chunk_size, IS_LENGTH - (size_t)tellg());
      } else {
        // If only asking for small amount of data, do not use full buffer
        // This is typically only called when header is read
        // Obviously will not trigger filling the next_file_buf
        load_from_file(max_bytes_to_decompress);
      }
    } else if(next_file_buf_cap > 0) {
      // Checks if need to swap file buffer
      swap_file_buffer_if_needed();
      // If above is triggered:
      // - Typically now file_buf_cap == FILE_BUFFER_CAP, file_buf_cursor == 0
      
      // NB file_buf_cap will only be < FILE_BUFFER_CAP if we have finished
      //   reading from the file, given next_file_buf is already filled with
      //   some data. In which case we would not be in this block.
      
      // Check if the file cursor has moved past the point where
      //   we should trigger filling more file chunks to next_file_buf
      
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
    if(file_buf_cap == file_buf_cursor) return(0);  // Finished reading file buffer
  }

  // No point asking for filling if next_file_buf is already at that level
  if(next_file_buf_cap >= spare_bytes_to_fill) spare_bytes_to_fill = 0;

  // Set decomp_threads = threads_to_use - 1 to use asynchronous file reading
  //   using the remaining thread
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
  
  bool check_gzip_head = false;
  uint32_t * u32; uint16_t * u16;

  // Trial run to see if data limited at source file cap, or destination cap
  unsigned int bgzf_count = 0;
  while(src_cursor < chunk_size) {
    if(file_buf_cursor + src_cursor + 28 > file_buf_cap) break;

    // Check valid bamGzipHead every 1000 blocks
    if(bgzf_count % 1000 == 0) check_gzip_head = true; // Check every 10000 blocks
    bgzf_count++;
    
    if(check_gzip_head) {
      if(strncmp(bamGzipHead, file_buf + file_buf_cursor + src_cursor, bamGzipHeadLength) != 0) {
        break;  // Gzip head corrupt
      } else {
        check_gzip_head = false;
      }
    }

    // Check entire bgzf block lies within file_buf
    u16 = (uint16_t*)(file_buf + file_buf_cursor + src_cursor + 16);
    if(file_buf_cursor + src_cursor + (*u16+1) > file_buf_cap) break;

    // Check entire bgzf block lies within chunk_size
    if(src_cursor + (*u16+1) > chunk_size) break;
    
    // Abort if filling dest with extra data will lead to overflow
    u32 = (uint32_t*)(file_buf + file_buf_cursor + src_cursor + (*u16+1) - 4);
    if(dest_cursor + *u32 > max_bytes_to_decompress) break;   
    
    // Increment cursors to next compressed / decompressed bgzf block
    src_cursor += *u16 + 1;
    dest_cursor += *u32;
  }

  if(check_gzip_head) {
    Rcpp::Rcout << "BGZF blocks corrupt\n";
    return(0);
  }

  // Stores how much data to decompress,
  //   and how much space to allocate to decompressed data
  size_t src_max = src_cursor; size_t dest_max = dest_cursor;
  src_cursor = 0; dest_cursor = 0;

// Divide the bgzf blocks among the number of decompression threads

  std::vector<size_t> src_bgzf_pos;   // source cursors
  std::vector<size_t> dest_bgzf_pos;  // dest cursors    
  std::vector<size_t> src_bgzf_cap;   // source cap: when thread should stop reading
  std::vector<size_t> dest_bgzf_cap;  // dest cap    

  // Ensures src_divider * decomp_threads > src_max
  size_t src_divider = 1 + (src_max / decomp_threads);
  size_t next_divider = std::min(src_divider, src_max);
  
  // Profile the file and increment until the destination buffer is reached
  unsigned int threads_accounted_for = 0;
  src_bgzf_pos.push_back(file_buf_cursor + src_cursor);  
  dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
  while(src_cursor < src_max) {
    // Increments cursors at next comp / decomp BGZF block
    u16 = (uint16_t*)(file_buf + file_buf_cursor + src_cursor + 16);
    u32 = (uint32_t*)(file_buf + file_buf_cursor + src_cursor + (*u16+1) - 4);
    src_cursor += *u16 + 1;
    dest_cursor += *u32;   

    // This block will trigger a maximum (decomp_threads - 1) times
    //   since src_divider * decomp_threads > src_max
    if(src_cursor > next_divider) {
      src_bgzf_cap.push_back(file_buf_cursor + src_cursor);  
      dest_bgzf_cap.push_back(decomp_cursor + dest_cursor);
      src_bgzf_pos.push_back(file_buf_cursor + src_cursor);  
      dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
      next_divider = std::min(next_divider + src_divider, src_max);
      threads_accounted_for++;
    }
  }
  // Now, src_cursor == src_max and dest_cursor == dest_max
  // Be pedantic and check this
  if(src_cursor != src_max || dest_cursor != dest_max) {
    Rcpp::Rcout << "Error occurred during BGZF block counting\n"
      << "This should not have occurred. Please report to mod author\n";
    return(0);
  }
  
  // threads_accounted_for == (decomp_threads - 1)
  //   but just in case it is not:
  while(threads_accounted_for < decomp_threads - 1) {
    src_bgzf_cap.push_back(file_buf_cursor + src_cursor);  
    dest_bgzf_cap.push_back(decomp_cursor + dest_cursor);
    src_bgzf_pos.push_back(file_buf_cursor + src_cursor);  
    dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
    threads_accounted_for++;
  }  

  // Fill last cap
  src_bgzf_cap.push_back(file_buf_cursor + src_max);   
  dest_bgzf_cap.push_back(decomp_cursor + dest_max);
   
  // Here, vector.size() == decomp_cursor
   
  // Now comes the multi-threaded decompression:
  bool error_occurred = false;
  
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads_to_use) schedule(static,1)
  #endif
  for(unsigned int k = 0; k < threads_to_use; k++) {
    if(k == decomp_threads) {
      // In asynchronous read, decomp_threads == threads_to_use - 1
      // Therefore, this last thread is used to prime spare file buffer
      read_file_chunk_to_spare_buffer(spare_bytes_to_fill);
      spare_bytes_to_fill = 0;
    } else {
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
            crc = crc32(crc32(0L, NULL, 0L), (Bytef*)(data_buf + thread_dest_cursor), *dest_size);
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
      }
    }
    
  }
  
  if(error_occurred) return(0);
  
  file_buf_cursor = src_bgzf_cap.at(src_bgzf_cap.size() - 1);
  data_buf_cap = dest_bgzf_cap.at(dest_bgzf_cap.size() - 1);

  return(dest_cursor);
}

// *************** Internal functions run by decompress() *********************

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

inline void pbam_in::move_residual_data(
    char** buf,                     // Pointer to char pointer (for file/data buffer)
    size_t & cursor, size_t & cap,  // Reference to buffer cursor / caps
    const size_t new_cap            // Desired size of new buffer
) {
  char * data_tmp;
  char * temp_buffer;
  size_t residual = cap - cursor;
  if(residual > 0) {
    // have to move bytes around
    temp_buffer = (char*)malloc(residual + 1);
    memcpy(temp_buffer, *buf + cursor, residual);
  
    *buf = (char*)realloc(data_tmp = *buf, 
      std::max(new_cap + 1, residual + 1));
    memcpy(*buf, temp_buffer, residual);
    free(temp_buffer);
    
    cap = residual;
  } else {
    *buf = (char*)realloc(data_tmp = *buf, new_cap + 1);
    cap = 0;
  }
  cursor = 0;
}

inline int pbam_in::clean_data_buffer(const size_t n_bytes_to_decompress) {
  move_residual_data(&data_buf, data_buf_cursor, data_buf_cap, n_bytes_to_decompress);
  return(0);
}

/* 
  If bytes remaining in primary file buffer is less than chunk_size,
  re-initialises file_buf, copying residual data from existing file_buf,
  and copies data from next_file_buf, upto FILE_BUFFER_CAP
*/
inline int pbam_in::swap_file_buffer_if_needed() {
  // Transfers residual data from file_buf to next_file_buf:
  if(next_file_buf_cap == 0) return(1);
  size_t chunk_size = (size_t)(FILE_BUFFER_CAP / chunks_per_file_buf);
  if(file_buf_cap  - file_buf_cursor > chunk_size) return(1);
  
  // Remove residual bytes to beginning of buffer:
  move_residual_data(&file_buf, file_buf_cursor, file_buf_cap, FILE_BUFFER_CAP);
  size_t residual = file_buf_cap;
  
  if(next_file_buf_cap <= FILE_BUFFER_CAP - file_buf_cap) {
    // If next_file_buf fits inside cap:
    //   transfer all contents over and free next_file_buf
    // Typically only happens when we have finished reading the BAM file

    memcpy(file_buf + file_buf_cap, next_file_buf, next_file_buf_cap);
    file_buf_cap += next_file_buf_cap;

    free(next_file_buf); next_file_buf = NULL;
    next_file_buf_cap = 0;
  } else {
    // Transfers some data from next_file_buf to fill up file_buf
    // Then move the residual in next_file_buf to beginning of buffer 
  
  // copies some data to fill file_buf to the brim
    memcpy(file_buf + file_buf_cap, next_file_buf, FILE_BUFFER_CAP - file_buf_cap);
    file_buf_cap = FILE_BUFFER_CAP;
    
    next_file_buf_cursor = FILE_BUFFER_CAP - residual;
    
    move_residual_data(&next_file_buf, next_file_buf_cursor, next_file_buf_cap, 
      FILE_BUFFER_CAP);
  }
  return(0);
}

/*  
  Read from file to fill n_bytes of file_buf
  First removes data upstream of file cursor
  Then fills buffer to n_bytes
*/
inline size_t pbam_in::load_from_file(const size_t n_bytes) {
  move_residual_data(&file_buf, file_buf_cursor, file_buf_cap, FILE_BUFFER_CAP);
  size_t residual = file_buf_cap;
  
  if(n_bytes <= residual) return(0);
  if(FILE_BUFFER_CAP <= file_buf_cap) return(0);

  size_t n_bytes_to_load =  std::min( std::max(n_bytes, residual) , FILE_BUFFER_CAP);  // Cap at file buffer
  size_t n_bytes_to_read = std::min(n_bytes_to_load - residual, IS_LENGTH - tellg());
  if(n_bytes_to_read == 0) return(0);
  
  read_file_to_buffer(file_buf + file_buf_cap, n_bytes_to_read);

  file_buf_cap += n_bytes_to_read;
  return(n_bytes_to_read);
}

inline size_t pbam_in::fill_file_buffer() {
  return(load_from_file(FILE_BUFFER_CAP));
}

/*  
  Read from file to fill n_bytes of next_file_buf
*/
inline size_t pbam_in::read_file_chunk_to_spare_buffer(const size_t n_bytes) {
  // There is no residual to move in next_file_buf
  size_t residual = next_file_buf_cap;

  if(n_bytes <= residual) return(0);
  if(FILE_BUFFER_CAP <= residual) return(0);
  
  size_t n_bytes_to_load =  std::min( std::max(n_bytes, residual) , FILE_BUFFER_CAP);  // Cap at file buffer
  size_t n_bytes_to_read = std::min(n_bytes_to_load - residual, IS_LENGTH - tellg());  
  if(n_bytes_to_read == 0) return(0);
  
  // Initialize next_file_buf here:
  if(next_file_buf_cap == 0) {
    char * file_tmp;
    next_file_buf = (char*)realloc(file_tmp = next_file_buf, FILE_BUFFER_CAP + 1);
  }
  
  read_file_to_buffer(next_file_buf + next_file_buf_cap, n_bytes_to_read);
  
  next_file_buf_cap += n_bytes_to_read;
  return(n_bytes_to_read);
}


#endif
