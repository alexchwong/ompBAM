#ifndef _ParaBAM_DEF
#define _ParaBAM_DEF

// [[Rcpp::depends(zlibbioc)]]
#include <zlib.h>
#include <zconf.h>

#include <fstream>    // std::ifstream

#include <cstring>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

union pb_header {
  char c[8];
  struct {
    char magic[4];
    int32_t l_text;
  } magic;
};

struct pb_core_32{
  // uint32_t block_size;
  int32_t refID;
  int32_t pos;
  uint8_t l_read_name;
  uint8_t mapq;
  uint16_t bin;
  uint16_t n_cigar_op;
  uint16_t flag;
  uint32_t l_seq;
  int32_t next_refID;
  int32_t next_pos;
  int32_t tlen;
};

class ParaBAM {
  private:
// ParaBAM settings
    size_t FILE_BUFFER_CAP = 100000000;     // 1e8
    size_t DATA_BUFFER_CAP = 1000000000;    // 1e9
    unsigned int FILE_BUFFER_SEGMENTS = 5;  // Divide file buffer into n segments
    unsigned int paraBAM_n_threads = 1;
  
    std::istream * IN;
    size_t IS_LENGTH;     // size of BAM file (Input Stream Length)
    
    size_t LAST_READ_BUFFER_POS;      // The position of the last incomplete read. 
    // There is no split read if LAST_READ_BUFFER_POS = BUFFER_SIZE
  
    // Header storage:
    char * magic_header;
    uint32_t l_text;
    char * headertext;
    uint32_t n_ref;
    std::vector<std::string> chr_names;
    std::vector<uint32_t> chr_lens;
  
    // Data buffers
    char * file_buf; size_t file_buf_cap; size_t file_buf_cursor;
    char * next_file_buf; size_t next_file_buf_cap;

    char * data_buf; size_t data_buf_cap; size_t data_buf_cursor;
    
    std::vector<char *> read_ptrs;
    std::vector<uint32_t> read_ptr_partitions;
    std::vector<uint32_t> read_cursors;

    size_t load_from_file(size_t n_bytes);
    size_t fill_file_buffer() {return(load_from_file(FILE_BUFFER_CAP));};
    size_t read_file_chunk_to_spare_buffer(size_t n_bytes);

    size_t decompress(size_t n_bytes_to_decompress);
    
    int swap_file_buffer_if_needed();
    unsigned int calculate_chunks_to_load_to_secondary_buffer();
    int clean_data_buffer(size_t n_bytes_to_decompress);

    unsigned int read(char * dest, unsigned int len);  // returns the number of bytes actually read
    unsigned int ignore(unsigned int len);
    unsigned int peek(char * dest, unsigned int len) ;

    bool fail() {return(IN->fail());};
    size_t tellg() {return((size_t)IN->tellg());};
  public:
    ParaBAM();
    ParaBAM(size_t file_buffer_cap, size_t data_buffer_cap, unsigned int file_buffer_segments);
    ~ParaBAM();

    int SetInputHandle(std::istream *in_stream, unsigned int n_threads_to_use); // opens file stream, checks EOF and file size
    int readHeader();
    int obtainChrs(std::vector<std::string> & s_chr_names, std::vector<uint32_t> & u32_chr_lens);

    size_t GetLength() { return(IS_LENGTH); };

    int fillReads();  // Returns 1 if no more reads available
    
    // Read functions:
    char * supplyRead(unsigned int thread_number = 0);
    pb_core_32 * readCore(char * read);
    char *readName(char * read, uint8_t & len);
    uint32_t *readCigar(char * read, uint16_t & len);
};

#include "ParaBAM_fns.h"

#endif

