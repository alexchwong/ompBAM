#ifndef _pbam_in_supplyRead
#define _pbam_in_supplyRead

inline pbam1_t pbam_in::supplyRead(const unsigned int thread_id) {
  pbam1_t read;
  if(thread_id > read_cursors.size()) {
    Rcpp::Rcout << "Invalid thread number parsed to supplyRead()\n";
    return(read);
  }
  if(read_cursors.at(thread_id) >= read_ptr_ends.at(thread_id)) {
    return(read);
  }
  read = pbam1_t(data_buf + read_cursors.at(thread_id), false);
  if(read.validate()) {
    read_cursors.at(thread_id) += read.block_size() + 4;
  } else {
    // Check this is actually end of thread buffer; throw error here otherwise
    if(read_cursors.at(thread_id) < read_ptr_ends.at(thread_id)) {
      Rcpp::Rcout << "Invalid read found before end of thread buffer " 
        << thread_id << ". read_cursor = " << read_cursors.at(thread_id)
        << ", read_ptr_ends = " << read_ptr_ends.at(thread_id) << '\n';
    }
  }    
  return(read);
}

#endif
