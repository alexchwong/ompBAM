#include "ParaBAM.h"

ParaBAM::~ParaBAM() {
  if(file_buf) free(file_buf);
  file_buf = NULL;
  if(data_buf) free(data_buf);
  data_buf = NULL;
}

int ParaBAM::openFile(std::istream *in_stream, unsigned int n_threads_to_use) {
  #ifdef _OPENMP
    if(n_threads_to_use > omp_get_max_threads()) {
      n_threads = omp_get_max_threads();
    } else {
      n_threads = n_threads_to_use;
    }
  #else
    n_threads = 1;
  #endif
  
  IN = in_stream;
  
  // Assign IS_LENGTH
  IN->seekg(0, std::ios_base::end);
  IS_LENGTH = IN->tellg();
  
  // Check valid EOF:
  IN->seekg(-bamEOFlength, std::ios_base::end);
  
  char eof_check[bamEOFlength + 1];
  IN->read(eof_check, bamEOFlength);
  if(strncmp(bamEOF, eof_check, bamEOFlength) != 0) {
    Rcout << "Error opening BAM - EOF bit corrupt. Perhaps this file is truncated?\n";
    IN = NULL;
    return(-1);
  }
  IN->clear(); IN->seekg(0, std::ios_base::beg);
  return(0);
}

size_t ParaBAM::load_from_file(size_t n_bytes_from_file) {
  /*  Loads file_buf from file
      Residual bytes from the previous load are transferred to the beginning
      of the buffer. Then the remaining capacity is the loaded.
      The total bytes available is stored in file_buf_cap.
      The read position is moved as appropriate.
      Returns: number of bytes actually read
  */
  char * temp_buffer;
  size_t residual = file_buf_cap-file_buf_cursor;
  // Remove residual bytes to beginning of buffer:
  if(residual > 0 && file_buf_cursor > 0) {
    // if file_buf_cursor == 0 then no move is required
    temp_buffer = (char*)malloc(residual + 1);
    memcpy(temp_buffer, file_buf + file_buf_cursor, residual);
  }
  
  char * file_tmp = file_buf;
  if(n_bytes_from_file <= residual) {
    // If requested bytes are less than residual data already in buffer
    // then clean and return the residual data. 0 bytes read from file
    file_buf = realloc(file_tmp, residual + 1);
    if(temp_buffer) memcpy(file_buf, temp_buffer, residual);
    free(temp_buffer);
    file_buf_cap = residual;
    file_buf_cursor = 0;
    return(0);
  }
  
  // Copies residual data to beginning of file
  file_buf = realloc(file_tmp, n_bytes_from_file + 1);
  if(temp_buffer) memcpy(file_buf, temp_buffer, residual);
  free(temp_buffer);
  file_buf_cursor = 0;
  
  // Read from file
  size_t = n_bytes_to_read = n_bytes_from_file - residual;
  size_t = n_bytes_remaining = IS_LENGTH - IN->tellg();
  if(n_bytes_to_read > n_bytes_remaining) n_bytes_to_read = n_bytes_remaining;
  IN->read(file_buf + residual, n_bytes_to_read);
  file_buf_cap = residual + n_bytes_to_read;
  return(n_bytes_to_read);
}

size_t ParaBAM::decompress(size_t n_bytes_to_decompress) {
/*
  Decompresses from the file buffer.
  Stops if the capacity 'n_bytes_to_decompress' is reached
  Also stops if all whole bgzf blocks from the file buffer are exhausted
  Returns the number of decompressed bytes added
*/
  
  // Remove residual bytes to beginning of buffer:
  char * temp_buffer;
  size_t residual = data_buf_cap-data_buf_cursor;
  
  if(residual > 0 && data_buf_cursor > 0) {
    // if data_buf_cursor == 0 then no move is required
    temp_buffer = (char*)malloc(residual + 1);
    memcpy(temp_buffer, data_buf + data_buf_cursor, residual);
  }
  
  char * data_tmp = data_buf;
  if(n_bytes_from_file <= residual) {
    // If requested bytes are less than residual data already in buffer
    // then clean and return the residual data. 0 bytes decompressed
    data_buf = realloc(data_tmp, residual + 1);
    if(temp_buffer) {
      memcpy(data_buf, temp_buffer, residual);
      free(temp_buffer);
    }    
    data_buf_cap = residual;
    data_buf_cursor = 0;
    return(0);
  }

  // Copies residual data to beginning of file
  data_buf = realloc(data_tmp, n_bytes_to_decompress + 1);
  if(temp_buffer) {
    memcpy(data_buf, temp_buffer, residual);
    free(temp_buffer);
  }
  data_buf_cursor = 0;

  size_t decomp_cursor = residual;  // The cursor to the data buffer to begin adding data
  size_t max_bytes_to_decompress = n_bytes_to_decompress - residual;
  // Check first gzip block:

  size_t src_cursor = 0;
  size_t dest_cursor = 0;
  unsigned int n_blocks_processed = 0;
  std::vector<size_t> src_bgzf_pos;
  std::vector<size_t> dest_bgzf_pos;
  bool check_gzip_head = true;
  uint32_t * u32; uint16_t * u16;

  // Profile the file and increment until the destination buffer is reached
  unsigned int bgzf_count = 0;
  while(1) {
    if(bgzf_count % 10000 == 1) check_gzip_head = true; // Check every 10000 blocks
    
    // Check bamGzipHead
    if(check_gzip_head) {
      if(strncmp(bamGzipHead, file_buf + file_buf_cursor, bamGzipHeadLength) != 0) {
        break;  // Gzip head corrupt
      } else {
        check_gzip_head = false;
      }
    }
    
    // Abort if cannot read smallest possible bgzf block:
    if(file_buf_cursor + src_cursor + 28 > file_buf_cap) break;
    u16 = (uint16_t*)(file_buf + file_buf_cursor + src_cursor + 16);
       
    // Abort if cannot read entire bgzf block:
    if(file_buf_cursor + src_cursor + (*u16+1) > file_buf_cap) break;
    u32 = (uint32_t*)(src + src_cursor + (*u16+1) - 4);
    
    // Abort if filling dest with extra data will lead to overflow
    if(dest_cursor + *u32 > max_bytes_to_decompress) break;

    // If checks passed, designate these blocks to be decompressed

    src_bgzf_pos.push_back(file_buf_cursor + src_cursor);   
    src_cursor += *u16 + 1;

    dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
    dest_cursor += *u32;

    bgzf_count++;
  }
  // NB: last src and dest bgzf pos are the end: So always one more bgzf pos than done
  src_bgzf_pos.push_back(file_buf_cursor + src_cursor);   
  dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
  
  // Assign jobs between n threads
  std::vector<unsigned int> thrd_partitions;
  unsigned int job_inc = 0;
  unsigned int job_size = (src_bgzf_pos.size()) / n_threads;
  for(unsigned int k = 0; k < n_threads; k++) {
    thrd_partitions.push_back(job_inc);
    job_inc+=job_size;
  }
  thrd_partitions.push_back(src_bgzf_pos.size()-1); 
  
  // Now comes the multi-threaded decompression:
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(unsigned int k = 0; k < n_threads; k++) {
  
    for(unsigned int i = thrd_partitions.at(k); i < thrd_partitions.at(k+1); i++) {
      size_t thread_src_cursor = 0;
      size_t thread_dest_cursor = 0;
      
      thread_src_cursor = src_bgzf_pos.at(i);
      thread_dest_cursor = dest_bgzf_pos.at(i);
      uint32_t * crc_check = (uint32_t *)(file_buf + src_bgzf_pos.at(i+1) - 8);
      uint32_t src_size = src_bgzf_pos.at(i+1) - thread_src_cursor;
      uint32_t dest_size = dest_bgzf_pos.at(i+1) - thread_dest_cursor;
      
      z_stream zs;
      zs.zalloc = NULL; zs.zfree = NULL; zs.msg = NULL;
      zs.next_in = (Bytef*)(file_buf + thread_src_cursor);
      zs.avail_in = src_size;
      zs.next_out = (Bytef*)(data_buf + thread_dest_cursor);
      zs.avail_out = dest_size;

      int ret = inflateInit2(&zs, -15);
      if(ret != Z_OK) {
          Rcout << "Exception during BAM decompression - inflateInit2() fail: (" << ret << ") "
            "BGZF block number" << i << '\n';
          return(ret);
      }
      ret = inflate(&zs, Z_FINISH);
      if(ret != Z_OK && ret != Z_STREAM_END) {
          std::ostringstream oss;
          Rcout << "Exception during BAM decompression - inflate() fail: (" << ret << ") "
            "BGZF block number" << i << '\n';
          return(ret);
      }
      ret = inflateEnd(&zs);
      // debug
      Rcout << "Bytes remaining to decompress = " << zs.avail_out << '\n';

      uint32_t crc = crc32(crc32(0L, NULL, 0L), 
        (Bytef*)(data_buf + thread_dest_cursor), dest_size);
      if(*crc_check != crc) {
          Rcout << "CRC fail during BAM decompression" << ", BGZF block number" << i << '\n';
          return(ret);
      }
    }
  
  }
  
  file_buf_cursor = src_bgzf_pos.at(src_bgzf_pos.size() - 1);
  data_buf_cursor = dest_bgzf_pos.at(dest_bgzf_pos.size() - 1);
  return(dest_cursor);
}

unsigned int ParaBAM::read(char * dest, unsigned int len) {
  // Reads bytes from the current data_buf, increments data_buff_cursor
  unsigned int n_bytes_to_read = min(len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}

unsigned int ParaBAM::ignore(unsigned int len) {
  // Ignores bytes from the current data_buf, increments data_buff_cursor
  unsigned int n_bytes_to_read = min(len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  // memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}

unsigned int ParaBAM::peek(char * dest, unsigned int len) {
  // Peeks at bytes from the current data_buf, does not increment data_buff_cursor
  unsigned int n_bytes_to_read = min(len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  // data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}

int ParaBAM::readHeader() {
  // Read 1 Mb from file:
  load_from_file(1000000);
  decompress(1000000);
  
  // Use conventional read to process header
  magic_header = (char*)malloc(8+1);
  read(magic_header, 8);
  if(strncpy(magic_header, magicstring, 4) != 0) {
    Rcout << "Invalid BAM magic string\n";
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
  
  char chrom_buffer[1000];
  for(unsigned int i = 0; i < n_ref; i++) {
    read(u32c, 4); u32 = (uint32_t *)u32c;
    read(chrom_buffer, *u32);
    chr_names.push_back(string(chrom_buffer, *u32-1));
    
    read(u32c, 4); u32 = (uint32_t *)u32c;
    chr_lens.push_back(*u32);
  }
  
  free(u32c);
  // We should be at end of reads.
}

int ParaBAM::obtainChrs(std::vector<std::string> & s_chr_names, std::vector<uint32_t> & u32_chr_lens) {
  if(n_ref == 0) {
    Rcout << "No chromosome names stored. Is ParaBAM::readHeader() been run yet?\n";
    return(0);
  }
  s_chr_names.clear();
  u32_chr_lens.clear();
  for(unsigned int i = 0; i < n_ref; i++) {
    s_chr_names.push_back(chr_names.at(i));
    u32_chr_lens.push_back(chr_lens.at(i));
  }
  return(n_ref);
}

