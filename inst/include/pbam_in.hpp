/* pbam_in.hpp pbam_in class

Copyright (C) 2021 Alex Chit Hei Wong

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.  */

#ifndef _pbam_in
#define _pbam_in

const size_t ompBAM_bufsize = 1024 * 1024 * 2;

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
      // Size of each of 2 file buffers for compressed data (default 500 Mb)
      const size_t data_buffer_cap, 
      // Size of data buffer for decompressed data (default 1 Gb)
      const unsigned int chunks_per_file_buffer, 
      // How many chunks per file buffer (default 5)
      const bool read_file_using_multiple_threads = true  
      // Whether to read threads in multi-threaded way (default true)
    );
    
    ~pbam_in();

    // File Openers //
    // NB: on file open, pbam_in automatically checks the BAM format and header

    // Opens a file directly for BAM reading. Supply a BAM file name as a string
    // Returns 0 if success, or -1 if error
    int openFile(std::string filename, unsigned int n_threads); 
    
    // Closes the file. Essentially clears the pbam_in buffers, except constructor settings
    int closeFile();
    
    /*
      Assign a ifstream handle to pbam_in; declares number of threads to use.
      in_stream must point to an ifstream that has opened a BAM file in binary mode
      Additionally, this function reads the BAM header.
    
      This is an alternate way to open a BAM file. The user opens the file
      themselves, using std::ios::in | std::ifstream::binary, then parses
      the ifstream handle to SetInputHandle

      Returns 0 if success, or -1 if error
    */
    int SetInputHandle(std::ifstream *in_stream, const unsigned int n_threads); 
   
    // Returns two vectors, containing chromosome names and lengths
    // Must be called after openFile() or SetInputHandle()
    // Returns the number of chromosomes, or -1 if error
    int obtainChrs(
      std::vector<std::string> & s_chr_names, 
      std::vector<uint32_t> & u32_chr_lens
    );

    /* 
      Reads the BAM file, decompressing to a maximum either by the data buffer cap,
        or the file chunk size (file_buf_cap / chunks_per_file_buf).
      Then assigns read cursors to point to the locations of reads to begin reading
        by each thread.
      If there were any reads that were not retrieved by supplyRead() since the
        last call to fillReads(), this function returns an error with message.
        
      Returns -1 if error, and 1 if EOF. Otherwise, returns 0
    */
    int fillReads();
    
    /* 
      Returns a pbam1_t of the next read to be returned by thread 
        (as specified by thread_id).
      This function is designed to be called within an OpenMP parallel loop, where:     
      - thread_id must be between [0, n_threads - 1].
      - supplyRead() should be called until it returns an invalid read.
      
      supplyRead() will return an invalid read if:
      - It is at end of buffer, or
      - The read buffer is corrupt
      
      The developer can distinguish between the two, by calling
        remainingThreadReadsBuffer(). This will return 0 if there are no
        threads remaining.
        
      Call pbam1_t::validate() to check whether the thread returned by
        this function is a "valid" read.
    */
    pbam1_t supplyRead(const unsigned int thread_id = 0);
    
    size_t remainingThreadReadsBuffer(const unsigned int thread_id = 0);
    
    // Returns the size of the opened BAM
    size_t GetFileSize() { return(IS_LENGTH); };

    // Returns the number of bytes decompressed
    size_t GetProgress() {return(prog_tellg());};
    
    int GetErrorState() {return(error_state);};
    /* 
      Returns the incremental number of bytes decompressed since the last call 
      to IncProgress(). A useful function for RcppProgress progress bars.
    */
    size_t IncProgress() {
      size_t INC = prog_tellg() - PROGRESS;
      PROGRESS = prog_tellg();
      return(INC);
    };

  private:
// pbam_in Settings
    size_t          FILE_BUFFER_CAP       = 5e8;
    size_t          DATA_BUFFER_CAP       = 1e9;
    unsigned int    chunks_per_file_buf   = 5;    // Divide file buffer into n segments
    unsigned int    threads_to_use        = 1;
    bool            multiFileRead         = true;
    std::string     FILENAME;
// File particulars
    std::ifstream    * IN;
    char INbuf[ompBAM_bufsize];
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
    size_t          next_file_buf_cursor;
    
    char *          data_buf; 
    size_t          data_buf_cap; 
    size_t          data_buf_cursor;
/* 
  Thread-specific read cursor positions and boundaries
  Thread returns a null read if cursor read_cursors >= read_ptr_ends
*/
    std::vector<size_t>         read_cursors;     // Cursor(s) of start of next read in each thread
    std::vector<size_t>         read_ptr_ends;    // Boundaries of read positions in each thread

// Error state of decompression
    int error_state = 0;


// Internal functions

    void            check_threads(unsigned int n_threads_to_check);
    int             check_file();

// *** Initialisers ***
    void            initialize_buffers();         // Initialises a pbam_in
    void            clear_buffers();              // Clears all buffers and re-initialises pbam_in

// *** Internal functions run by decompress() ***

    int             read_file_to_buffer(char * buf, const size_t len);

    // Moves residual data to beginning of file or data buffer
    // Buffer is passed by reference to the pointer
    void move_residual_data(
        char** buf,                     // Pointer to char pointer (for file/data buffer)
        size_t & cursor, size_t & cap,  // Reference to buffer cursor / caps
        const size_t new_cap            // Desired size of new buffer
    );

    // Removes all data upstream of data_buf_cursor. Then, sets data_buf_cursor to 0.
    int             clean_data_buffer(const size_t n_bytes_to_decompress);  
    
    /* 
      If bytes remaining in primary file buffer is less than chunk_size,
      re-initialises file_buf, copying residual data from existing file_buf,
      and copies data from next_file_buf, upto FILE_BUFFER_CAP
    */
    int             swap_file_buffer_if_needed();
    
    // Reads from file and copies data to file_buf, up to a maximum of n_bytes
    size_t          load_from_file(const size_t n_bytes);

    // Runs load_from_file using n_bytes = FILE_BUFFER_CAP
    size_t          fill_file_buffer();

    // Essentially same as load_from_file, but this reads data into next_file_buf
    size_t          read_file_chunk_to_spare_buffer(const size_t n_bytes);

    /* 
      Main Decompress Function
      - Returns the number of bytes successfully decompressed
    */
    size_t          decompress(const size_t n_bytes_to_decompress);

// *** Internal functions used by readHeader() ***
    unsigned int read(char * dest, const unsigned int len);  // returns the number of bytes actually read
    unsigned int ignore(unsigned const int len);
    unsigned int peek(char * dest, const unsigned int len) ;
    
// *** Reads the BAM header. Automatically run with SetInputHandle() or openFile() ***
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

#include "pbam_in_constructors.hpp"
#include "pbam_in_IO.hpp"
#include "pbam_in_decompress.hpp"
#include "pbam_in_fillReads.hpp"
#include "pbam_in_supplyRead.hpp"
#include "pbam_in_internals.hpp"

#endif