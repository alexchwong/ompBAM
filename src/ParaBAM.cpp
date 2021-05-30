#include "ParaBAM.h"
#include "Rcpp.h"
using namespace Rcpp;

static const int bamGzipHeadLength = 16;  // +2 a uint16 with the full block length.
static const char bamGzipHead[bamGzipHeadLength+1] = 
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00";
static const int bamEOFlength = 28;
static const char bamEOF[bamEOFlength+1] =
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";

static const int magiclength = 4;
static const char magicstring[magiclength+1] = "\x42\x41\x4d\x01";

ParaBAM::ParaBAM() {
  file_buf = NULL;
  data_buf = NULL;
  file_buf_cap = 0; file_buf_cursor = 0;
  data_buf_cap = 0; data_buf_cursor = 0;
  next_file_buf_cap = 0;
  paraBAM_n_threads = 1;
  magic_header = NULL; l_text = 0; headertext = NULL; n_ref = 0;
  
  read_ptrs.resize(0);
  read_ptr_partitions.resize(0);
  read_cursors.resize(0);
  
  magic_header = NULL;
  l_text = 0;
  headertext = NULL;
  n_ref = 0;
  chr_names.resize(0);
  chr_lens.resize(0);
  
}

ParaBAM::ParaBAM(size_t file_buffer_cap, size_t data_buffer_cap, unsigned int file_buffer_segments) {
  file_buf = NULL;
  data_buf = NULL;
  file_buf_cap = 0; file_buf_cursor = 0;
  data_buf_cap = 0; data_buf_cursor = 0;
  next_file_buf_cap = 0;
  paraBAM_n_threads = 1;
  magic_header = NULL; l_text = 0; headertext = NULL; n_ref = 0;
  
  read_ptrs.resize(0);
  read_ptr_partitions.resize(0);
  read_cursors.resize(0);
  
  magic_header = NULL;
  l_text = 0;
  headertext = NULL;
  n_ref = 0;
  chr_names.resize(0);
  chr_lens.resize(0);
  
  FILE_BUFFER_CAP = file_buffer_cap;
  DATA_BUFFER_CAP = data_buffer_cap;
  FILE_BUFFER_SEGMENTS = file_buffer_segments;
  // TODO: sanity checks for these parameters
}

ParaBAM::~ParaBAM() {
  if(file_buf) free(file_buf);
  file_buf = NULL;
  if(data_buf) free(data_buf);
  data_buf = NULL;
  if(magic_header) free(magic_header);
  magic_header = NULL;
  if(headertext) free(headertext);
  headertext = NULL;
}

int ParaBAM::SetInputHandle(std::istream *in_stream, unsigned int n_threads_to_use) {
  #ifdef _OPENMP
    if(n_threads_to_use > (unsigned int)omp_get_max_threads()) {
      paraBAM_n_threads = (unsigned int)omp_get_max_threads();
    } else {
      paraBAM_n_threads = n_threads_to_use;
    }
  #else
    paraBAM_n_threads = 1;
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
  return(0);
}

int ParaBAM::swap_file_buffer() {
  // Transfers residual data from current file buffer to the next:
  
  char * file_tmp;
  char * residual_data_buffer;
  size_t residual = file_buf_cap-file_buf_cursor;
  // Remove residual bytes to beginning of buffer:
  if(residual > 0) {
    // have to move bytes around
    residual_data_buffer = (char*)malloc(residual + 1);
    memcpy(residual_data_buffer, file_buf + file_buf_cursor, residual);
  
    file_buf = (char*)realloc(file_tmp = file_buf, FILE_BUFFER_CAP);
    memcpy(file_buf, residual_data_buffer, residual);
    free(residual_data_buffer);

  } else {
    // residual == 0
    file_buf = (char*)realloc(file_tmp = file_buf, FILE_BUFFER_CAP);
  }

  // Transfers any daata from next_file_buf to file_buf:
  if(residual + next_file_buf_cap > FILE_BUFFER_CAP) {
    file_buf = (char*)realloc(file_tmp = file_buf, residual + next_file_buf_cap + 1);
  }
  memcpy(file_buf + residual, next_file_buf, next_file_buf_cap);

  // Reset file stuff:
  file_buf_cap = next_file_buf_cap + residual;
  file_buf_cursor = 0;

  file_tmp = next_file_buf;
  next_file_buf = (char*)realloc(file_tmp = next_file_buf, FILE_BUFFER_CAP);
  next_file_buf_cap = 0;
  
  return(0);
}


size_t ParaBAM::load_from_file(size_t n_bytes) {
  // Read from file to fill n_bytes of file_buf
  // Rcout << "load_from_file()\n";
  if(n_bytes < file_buf_cap) return(0);
  size_t n_bytes_to_read = n_bytes - file_buf_cap;
  size_t n_bytes_remaining = IS_LENGTH - IN->tellg();

  if(n_bytes_to_read > n_bytes_remaining) n_bytes_to_read = n_bytes_remaining;
  
  char * file_tmp = file_buf;
  file_buf = (char*)realloc(file_tmp, file_buf_cap + n_bytes_to_read + 1);
  
  IN->read(file_buf + file_buf_cap, n_bytes_to_read);
  file_buf_cap += n_bytes_to_read;
  // Rcout << "load_from_file() done\n";
  return(n_bytes_to_read);
}

size_t ParaBAM::read_file_chunk_to_spare_buffer(size_t n_bytes) {
  // Read from file to fill n_bytes of next_file_buf
  if(n_bytes < next_file_buf_cap) return(0);
  size_t n_bytes_to_read = n_bytes - next_file_buf_cap;
  size_t n_bytes_remaining = IS_LENGTH - IN->tellg();

  if(n_bytes_to_read > n_bytes_remaining) n_bytes_to_read = n_bytes_remaining;
  
  if(next_file_buf_cap + n_bytes_to_read > FILE_BUFFER_CAP) {
    char * file_tmp;
    next_file_buf = (char*)realloc(file_tmp = next_file_buf, file_buf_cursor + n_bytes_to_read + 1);
  }
  IN->read(next_file_buf + next_file_buf_cap, n_bytes_to_read);
  next_file_buf_cap += n_bytes_to_read;
  return(n_bytes_to_read);
}

int ParaBAM::clean_data_buffer(size_t n_bytes_to_decompress) {
  // Remove residual bytes to beginning of buffer:
  // Rcout << "clean_data_buffer\n";
  char * data_tmp = data_buf;
  char * temp_buffer;
  size_t residual = data_buf_cap-data_buf_cursor;
  
  if(residual > 0) {
    // have to move bytes around

    temp_buffer = (char*)malloc(residual + 1);
    memcpy(temp_buffer, data_buf + data_buf_cursor, residual);
  
    data_buf = (char*)realloc(data_tmp, max(n_bytes_to_decompress + 1, residual + 1));
    memcpy(data_buf, temp_buffer, residual);
    free(temp_buffer);
    
    data_buf_cap = residual;

  } else {
    // residual == 0
    data_buf = (char*)realloc(data_tmp, n_bytes_to_decompress + 1);
    data_buf_cap = 0;
  }
  data_buf_cursor = 0;
  // Rcout << "clean_data_buffer done\n";

  return(0);
}

size_t ParaBAM::decompress(size_t n_bytes_to_decompress) {
/*
  Wipes all data behind data_buf_cursor
  If file_buf data remaining is below residual (FILE_BUFFER_CAP / FILE_BUFFER_SEGMENTS)
    then runs swap_file_buffer
  If file_buf_cursor / FILE_BUFFER_SEGMENTS > next_file_buf_cap / FILE_BUFFER_SEGMENTS
    then runs read_file_chunk_to_spare_buffer
  Decompresses any data in file_buf to fill up to n_bytes_to_decompress in data_buf
*/

  // Rcout << "decompress()\n";

  clean_data_buffer(n_bytes_to_decompress);

  // Rcout << "data_buf_cap " << data_buf_cap << " n_bytes_to_decompress " << n_bytes_to_decompress << "\n";
  if(data_buf_cap >= n_bytes_to_decompress) return(0);

// Check if file buffer needs filling
  size_t spare_bytes_to_fill = 0;
  size_t chunk_size = (size_t)(FILE_BUFFER_CAP / FILE_BUFFER_SEGMENTS);
  if(!IN->eof()) {
    // If asking to decompress a large amount of data:
    if(n_bytes_to_decompress > chunk_size) {

      if(next_file_buf_cap == 0 && file_buf_cap < 
          FILE_BUFFER_CAP - chunk_size) {
        // If primary buffer is underfilled, do this first:
        fill_file_buffer();
      } else if(next_file_buf_cap > 0 && file_buf_cap - file_buf_cursor < chunk_size) {
        // Checks if need to swap file buffer 
        swap_file_buffer();
      }
      unsigned int n_chunks_to_fetch = (file_buf_cursor / chunk_size) - 
          (next_file_buf_cap - chunk_size);
      if(n_chunks_to_fetch > 0) {
        spare_bytes_to_fill = chunk_size * (unsigned int)(file_buf_cursor / chunk_size) +
          chunk_size;
      }
    } else if(n_bytes_to_decompress > (file_buf_cap - file_buf_cursor)) {
      load_from_file(n_bytes_to_decompress + file_buf_cursor);
    }
    // TODO: account for if user asks repeatedly for small chunks of data to a point
    //   where it overfills the file buffer
  }

  size_t decomp_cursor = data_buf_cap;  // The cursor to the data buffer to begin adding data
  size_t max_bytes_to_decompress = n_bytes_to_decompress - decomp_cursor;
  // Rcout << "max_bytes_to_decompress " << max_bytes_to_decompress << '\n';
  size_t src_cursor = 0;
  size_t dest_cursor = 0;
  
  std::vector<size_t> src_bgzf_pos;
  std::vector<size_t> dest_bgzf_pos;
  bool check_gzip_head = true;
  uint32_t * u32; uint16_t * u16;

  // Rcout << "file_buf_cursor = " << file_buf_cursor << '\n';

  // Profile the file and increment until the destination buffer is reached
  unsigned int bgzf_count = 0;
  while(1) {
    if(bgzf_count % 10000 == 1) check_gzip_head = true; // Check every 10000 blocks
    
    // Check bamGzipHead
    if(check_gzip_head) {
      if(strncmp(bamGzipHead, file_buf + file_buf_cursor + src_cursor, bamGzipHeadLength) != 0) {
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
    u32 = (uint32_t*)(file_buf + file_buf_cursor + src_cursor + (*u16+1) - 4);
    // Rcout << "*u32 " << *u32 << '\n';
    // Abort if filling dest with extra data will lead to overflow
    if(dest_cursor + *u32 > max_bytes_to_decompress) break;

    // If checks passed, designate these blocks to be decompressed

    src_bgzf_pos.push_back(file_buf_cursor + src_cursor);  
    src_cursor += *u16 + 1;
    
    dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
    dest_cursor += *u32;

    bgzf_count++;
  }
  // Rcout << "dest_cursor " << dest_cursor << '\n';
  if(check_gzip_head) {
    Rcout << "Gzip header corrupt\n";
    return(0);
  }
  
  // NB: last src and dest bgzf pos are the end: So always one more bgzf pos than done
  src_bgzf_pos.push_back(file_buf_cursor + src_cursor);   
  dest_bgzf_pos.push_back(decomp_cursor + dest_cursor);
  
  unsigned int decomp_threads = paraBAM_n_threads;
  if(spare_bytes_to_fill > 0) {
    if(decomp_threads > 1) decomp_threads--;
  } else {
    read_file_chunk_to_spare_buffer(spare_bytes_to_fill);
    spare_bytes_to_fill = 0;
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
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(unsigned int k = 0; k < paraBAM_n_threads; k++) {
    if(k == decomp_threads) {
      // This only occurs if spare_bytes_to_fill > 0
      
      // Last thread is used to prime spare file buffer
      read_file_chunk_to_spare_buffer(spare_bytes_to_fill);
      spare_bytes_to_fill = 0;
    } else {

      unsigned int i = thrd_partitions.at(k);
      bool error_occurred = false;
      while(i < thrd_partitions.at(k+1) && !error_occurred) {
        size_t thread_src_cursor = 0;
        size_t thread_dest_cursor = 0;
        
        thread_src_cursor = src_bgzf_pos.at(i);
        thread_dest_cursor = dest_bgzf_pos.at(i);
        
        uint32_t * crc_check = (uint32_t *)(file_buf + src_bgzf_pos.at(i+1) - 8);
        uint32_t src_size = src_bgzf_pos.at(i+1) - thread_src_cursor;
        uint32_t dest_size = dest_bgzf_pos.at(i+1) - thread_dest_cursor;
        // Rcout << thread_src_cursor << "\t" << thread_dest_cursor << '\t'
          // << src_size << "\t" << dest_size << '\n';
        z_stream zs;
        zs.zalloc = NULL; zs.zfree = NULL; zs.msg = NULL;
        zs.next_in = (Bytef*)(file_buf + thread_src_cursor + 18);
        zs.avail_in = src_size - 18;
        zs.next_out = (Bytef*)(data_buf + thread_dest_cursor);
        zs.avail_out = dest_size;

        int ret = inflateInit2(&zs, -15);
        if(ret != Z_OK) {
            Rcout << "Exception during BAM decompression - inflateInit2() fail: (" << ret << ") "
              "BGZF block number" << i << '\n';
            error_occurred = true;
        }
        if(!error_occurred) {
          ret = inflate(&zs, Z_FINISH);
          if(ret != Z_OK && ret != Z_STREAM_END) {
              Rcout << "Exception during BAM decompression - inflate() fail: (" << ret << ") "
                "BGZF block number" << i << '\n';
              error_occurred = true;
          }
          ret = inflateEnd(&zs);
          // debug
          // Rcout << "Bytes remaining to decompress = " << zs.avail_out << '\n';

          if(!error_occurred) {
            uint32_t crc = crc32(crc32(0L, NULL, 0L), 
              (Bytef*)(data_buf + thread_dest_cursor), dest_size);
            if(*crc_check != crc) {
                Rcout << "CRC fail during BAM decompression" << ", BGZF block number" << i << '\n';
                error_occurred = true;
            }
          }
        }

        i++;
      }
    }
    
  }
  
  file_buf_cursor = src_bgzf_pos.at(src_bgzf_pos.size() - 1);
  data_buf_cap = dest_bgzf_pos.at(dest_bgzf_pos.size() - 1);
  // Rcout << "file_buf_cursor = " << file_buf_cursor << '\n';
  // Rcout << "decompress() done\n";
  return(dest_cursor);
}

unsigned int ParaBAM::read(char * dest, unsigned int len) {
  // Reads bytes from the current data_buf, increments data_buff_cursor
  unsigned int n_bytes_to_read = min((size_t)len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  // Rcout << n_bytes_to_read << '\n';
  memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}

unsigned int ParaBAM::ignore(unsigned int len) {
  // Ignores bytes from the current data_buf, increments data_buff_cursor
  unsigned int n_bytes_to_read = min((size_t)len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  // memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}

unsigned int ParaBAM::peek(char * dest, unsigned int len) {
  // Peeks at bytes from the current data_buf, does not increment data_buff_cursor
  unsigned int n_bytes_to_read = min((size_t)len, data_buf_cap - data_buf_cursor);
  if(n_bytes_to_read == 0) return(0);
  memcpy(dest, data_buf + data_buf_cursor, n_bytes_to_read);
  // data_buf_cursor += n_bytes_to_read;
  return(n_bytes_to_read); 
}

int ParaBAM::readHeader() {
  // Read 1 Mb from file:
  decompress(1000000);

  magic_header = (char*)malloc(8+1);
  read(magic_header, 8);
  if(strncmp(magic_header, magicstring, 4) != 0) {
    Rcout << "Invalid BAM magic string\n";
    return(-1);
  }
  uint32_t * u32 = (uint32_t *)(magic_header + 4);
  l_text = *u32;
  if(l_text > 1000000) {
    decompress(l_text);
  }
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
    chrName = string(chrom_buffer, *u32-1);
    chr_names.push_back(chrName);
    // Rcout << " chr_names " << chrName << '\n';
    read(u32c, 4); u32 = (uint32_t *)u32c;
    chr_lens.push_back(*u32);
  }
  
  free(u32c);
  // We should be at end of reads.
  return(0);
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

int ParaBAM::fillReads() {
  // Returns -1 if error, and 1 if EOF. Otherwise, returns 0
  Rcout << "fillReads()\n";
  // Returns -1 if header not read:
  if(chr_names.size() == 0) {
    Rcout << "Header has not been read\n";
    return(-1);
  } 
  
  // Check if previous reads are all read:
  if(read_cursors.size() > 0) {
    for(unsigned int i = 0; i < read_cursors.size(); i++) {
      if(read_cursors.at(i) < read_ptr_partitions.at(i)) {
        Rcout << "Thread " << i << " has reads remaining. Please debug your code "
          << "and make sure all threads clear their reads before filling any more reads\n";
        return(-1);
      }
    }
  }
  
  // Clear read pointers:
  read_ptrs.resize(0);
  read_ptr_partitions.resize(0);
  read_cursors.resize(0);

  decompress(DATA_BUFFER_CAP);
  // if(ret == 0) {
    // Rcout << "End of buffer reached\n";
    // return(1);
  // }
  
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
    Rcout << "End of buffer reached\n";
    return(1);
  }
  
  // Partition reads by number of threads
  unsigned int reads_per_thread = 1 + (read_ptrs.size() / paraBAM_n_threads);
  unsigned int cursor = 0;
  while(read_cursors.size() < paraBAM_n_threads) {
    read_cursors.push_back(cursor);
    cursor += reads_per_thread;
    if(cursor > read_ptrs.size()) cursor = read_ptrs.size();
    read_ptr_partitions.push_back(cursor);
  }
  return(0);
}

char * ParaBAM::supplyRead(unsigned int thread_number) {
  if(thread_number > read_cursors.size()) {
    Rcout << "Invalid thread number parsed to supplyRead()\n";
    return(NULL);
  }
  if(read_cursors.at(thread_number) == read_ptr_partitions.at(thread_number)) {
    return(NULL);
  }
  char * read = read_ptrs.at(read_cursors.at(thread_number));
  read_cursors.at(thread_number)++;
  
  // pb_core_32 * core = (pb_core_32 *)(read + 4);
  // Rcout << "Read refID = " << core->refID << ", pos = " << core->pos << '\n';
  return(read);
}

pb_core_32 * ParaBAM::readCore(char * read) {
  pb_core_32 *core = (pb_core_32*)(read + 4);
  return(core);
}

char * ParaBAM::readName(char * read, uint8_t & len) {
  if(!read) return(NULL);
  pb_core_32 * core = (pb_core_32 *)(read + 4);
  
  len = core->l_read_name;
  char *read_name = (char*)(read + 36);
  return(read_name);
}

uint32_t * ParaBAM::readCigar(char * read, uint16_t & len) {
  if(!read) return(NULL);
  pb_core_32 * core = (pb_core_32 *)(read + 4);
  
  len = core->n_cigar_op;
  uint32_t * cigar = (uint32_t*)(read + 36 + core->l_read_name);
  return(cigar);
}