#ifndef _pbam_in_supplyRead
#define _pbam_in_supplyRead

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
