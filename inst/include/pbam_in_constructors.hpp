/* pbam_in_constructors.hpp pbam_in class constructors

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

#ifndef _pbam_in_constructors
#define _pbam_in_constructors

// Declare pbam_in with default settings
inline pbam_in::pbam_in() {
  initialize_buffers();
}

// Declare pbam_in with custom settings; for advanced users only
inline pbam_in::pbam_in(
  const size_t file_buffer_cap,   // Default 500 Mb
  const size_t data_buffer_cap,   // Default 1 Gb
  const unsigned int chunks_per_file_buffer,    // Default 5 - i.e. 40 Mb file chunks
  // File read triggers after each 40 Mb decompressed
  const bool read_file_using_multiple_threads   // default true
) {
  initialize_buffers();
  
  // Limits chunk size to be above 1 Mb
  if(file_buffer_cap / chunks_per_file_buffer < 1024576) {
    cout << "FILE_BUFFER_CAP / chunks_per_file_buffer (chunk size) must be above 1Mb\n";
    return;
  }
  // Ensures decompressed data buffer is bigger than file buffer. For sanity reasons.
  if(data_buffer_cap < file_buffer_cap) {
    cout << "DATA_BUFFER_CAP must not be smaller than FILE_BUFFER_CAP\n";
    return;
  }  

  FILE_BUFFER_CAP = file_buffer_cap;
  DATA_BUFFER_CAP = data_buffer_cap;
  chunks_per_file_buf = chunks_per_file_buffer;
  threads_to_use = 1;
  multiFileRead = read_file_using_multiple_threads;
}

inline pbam_in::~pbam_in() {
  clear_buffers();
}

#endif
