/* pbam_in_supplyRead.hpp pbam_in supplyRead()

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

#ifndef _pbam_in_supplyRead
#define _pbam_in_supplyRead

inline pbam1_t pbam_in::supplyRead(const unsigned int thread_id) {
  pbam1_t read;
  if(thread_id > read_cursors.size()) {
    cout << "Invalid thread number parsed to supplyRead()\n";
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
      cout << "Invalid read found before end of thread buffer " 
        << thread_id << ". read_cursor = " << read_cursors.at(thread_id)
        << ", read_ptr_ends = " << read_ptr_ends.at(thread_id) << '\n';
    }
  }    
  return(read);
}

#endif
