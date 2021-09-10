#ifndef _pbam_in_fillReads
#define _pbam_in_fillReads


inline int pbam_in::fillReads() {
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
        Rcpp::Rcout << "Thread " << i << " has reads remaining. Please debug your code "
          << "and make sure all threads clear their reads before filling any more reads\n";
        return(-1);
      }
    }
  }
  
  // Clear read pointers:
  read_cursors.resize(0);
  read_ptr_ends.resize(0);
  
  // Call decompress
  decompress(DATA_BUFFER_CAP);
  
  // Check decompressed data contains at least 1 full read  
  if(data_buf_cap - data_buf_cursor < 4) return(1);
  uint32_t *u32p;
  u32p = (uint32_t *)(data_buf + data_buf_cursor);
  if(*u32p + 4 > data_buf_cap - data_buf_cursor) return(1);

  // Roughly divide the buffer into N regions:
  size_t data_divider = 1 + ((data_buf_cap - data_buf_cursor) / threads_to_use);
  size_t next_divider = std::min(data_buf_cursor + data_divider, data_buf_cap);

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
    if(data_buf_cursor >= next_divider) {
      read_ptr_ends.push_back(data_buf_cursor);
      read_cursors.push_back(data_buf_cursor);
      next_divider = std::min(data_buf_cap, next_divider + data_divider);
      threads_accounted_for++;
    }
  }

  while(threads_accounted_for < threads_to_use - 1) {
    read_ptr_ends.push_back(data_buf_cursor);
    read_cursors.push_back(data_buf_cursor);
    threads_accounted_for++;
  }
  read_ptr_ends.push_back(data_buf_cursor);

  return(0);
}

// Internal
inline size_t pbam_in::remainingThreadReadsBuffer(const unsigned int thread_number) {
  if(thread_number > threads_to_use) {
    Rcpp::Rcout << "pbam_in object was not initialized with " << 
      thread_number << " threads\n";
    return(0);
  }
  if(read_cursors.size() <= thread_number) {
    Rcpp::Rcout << "Thread " << thread_number << " is not initialized with reads\n";
    return(0);
  }
  return(read_ptr_ends.at(thread_number) - read_cursors.at(thread_number));
}

#endif
