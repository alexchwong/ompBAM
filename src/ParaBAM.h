#ifndef _ParaBAM
#define ParaBAM

#include <zlib.h>
// #include <zconf.h>
#include <iostream>
#include <fstream>    // std::ifstream
#include <algorithm>  // std::sort
#include <functional> // std::function

#include <cstring>
#include <vector>
#include <map>
#include <sys/types.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef Rcout
#define Rcout cout
#endif

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

static const int bamGzipHeadLength = 16;  // +2 a uint16 with the full block length.
static const char bamGzipHead[bamGzipHeadLength+1] = 
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00";
static const int bamEOFlength = 28;
static const char bamEOF[bamEOFlength+1] =
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";

static const int magiclength = 4;
static const char magicstring[magiclength+1] = "\x42\x41\x4d\x01";


class ParaBAM {
  private:
    unsigned int n_threads = 1;
  
    istream * IN;
    size_t IS_LENGTH;     // size of BAM file (Input Stream Length)
    
    size_t LAST_READ_BUFFER_POS = 0;      // The position of the last incomplete read. 
    // There is no split read if LAST_READ_BUFFER_POS = BUFFER_SIZE
  
    // Header storage:
    char * magic_header;
    uint32_t l_text = 0;
    char * headertext;
    uint32_t n_ref = 0;
    std::vector<string> chr_names;
    std::vector<uint32_t> chr_lens;
  
    char * file_buf; size_t file_buf_cap = 0; size_t file_buf_cursor = 0;
    char * data_buf; size_t data_buf_cap = 0; size_t data_buf_cursor = 0
    
    std::vector<pb_core_32 *> read_ptrs;
    std::vector<uint32_t> read_cursors;
    
    size_t load_from_file(size_t n_bytes_from_file);
    size_t decompress(size_t n_bytes_to_decompress);
    
  public:
    ~ParaBAM();

    int openFile(std::istream *in_stream, unsigned int n_threads_to_use); // opens file stream, checks EOF and file size
    int readHeader();
    int obtainChrs(std::vector<std::string> & s_chr_names, std::vector<uint32_t> & u32_chr_lens);


    unsigned int read(char * dest, unsigned int len);  // returns the number of bytes actually read
    unsigned int ignore(unsigned int len);
    unsigned int peek(char * dest, unsigned int len) ;

    bool fail() {return(IN->fail());};
    size_t tellg() {return((size_t)IN->tellg());};
    size_t GetLength() { return(IS_LENGTH); };
    
    
}


#endif
