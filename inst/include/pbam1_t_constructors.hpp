/* pbam1_t_constructors.hpp pbam1_t constructors

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

#ifndef _pbam1_t_constructors
#define _pbam1_t_constructors

// ************************** Constructor (empty) ******************************
inline pbam1_t::pbam1_t() {
  read_buffer = NULL;
  realized = false;
  core = NULL;
  block_size_val = 0;   tag_size_val = 0;
}

// ********************************* Destructor ********************************
inline pbam1_t::~pbam1_t() {
  if(read_buffer && realized) {
    free(read_buffer);
    read_buffer = NULL;
  }
  realized = false;
  core = NULL;
  block_size_val = 0;   tag_size_val = 0;
}

// ********************************* Validity **********************************

inline bool pbam1_t::validate() const {
  if(!read_buffer) return(false);

  uint32_t *temp_block_size = (uint32_t *)(read_buffer);
  if(*temp_block_size != block_size_val) return(false);
  if(!core) return(false);
  
  // Quick hash to check core; assume buffer is valid if core is valid
  if(tag_size_val != block_size_val - (
    32 + core->l_read_name + core->n_cigar_op * 4 + 
        core->l_seq + ((core->l_seq + 1) / 2)
  )) {
    // Attempt to identify which read is causing the fail:
    std::string read_name_s;
    char *tmp = (char*)(read_buffer + 36);
    read_name_s.assign(tmp);
    cout << "Invalid read: " << read_name_s << "\n";
    return(false);
  }
  return(true);
}

// ************************** Copy read to real buffer *************************

inline int pbam1_t::realize() {
  if(realized) return(0);
  if(validate()) {
    char *tmp = read_buffer;
    read_buffer = (char*)malloc(block_size_val + 5);
    memcpy(read_buffer, tmp, block_size_val + 4);
    core = (pbam_core_32 *)(read_buffer + 4);
    
    uint32_t *temp_block_size = (uint32_t *)(read_buffer);
    block_size_val = *temp_block_size;
    // re-derive block_size_val and tag_size_val:
    tag_size_val = block_size_val - (32 + 
        core->l_read_name + core->n_cigar_op * 4 + 
        core->l_seq + ((core->l_seq + 1) / 2));
    
    realized = true;
  }
  if(!validate()) return(-1);
  return(0);
}

// ******** Constructor given pointer to buffer; option to realize read ********
// This function is called by pbam_in::supplyRead()
inline pbam1_t::pbam1_t(char * src, bool realize) {
  
  
  uint32_t *temp_block_size = (uint32_t *)(src);
  block_size_val = *temp_block_size;
  read_buffer = src;
  // temp assignment to direct buffer
  core = (pbam_core_32 *)(read_buffer + 4);
  
  // Check tag size is zero or positive_sign
  if(32 + core->l_read_name + core->n_cigar_op * 4 + 
        core->l_seq + ((core->l_seq + 1) / 2) > block_size_val) {
    reset();
  } else {
    tag_size_val = block_size_val - (32 + 
        core->l_read_name + core->n_cigar_op * 4 + 
        core->l_seq + ((core->l_seq + 1) / 2));
    if(realize) {
      read_buffer = (char*)malloc(block_size_val + 1);
      memcpy(read_buffer, src, block_size_val);      
      realized = true;
    } else {
      read_buffer = src;
      realized = false;
    }
    core = (pbam_core_32 *)(read_buffer + 4);
    validate();
  }
}

// ********************************* Internals *********************************

// *********************************** Reset ***********************************
inline void pbam1_t::reset() {
  if(read_buffer && realized) {
    free(read_buffer);
    read_buffer = NULL;
  }
  realized = false;
  core = NULL;
  block_size_val = 0;   tag_size_val = 0;
}


// ***************************** Copy constructor *****************************
inline pbam1_t::pbam1_t(const pbam1_t &t) {
  if(t.isReal()) {
    if(t.validate()) {
      read_buffer = (char*)malloc(t.block_size_val + 1);
      memcpy(read_buffer, t.read_buffer, t.block_size_val);
      block_size_val = t.block_size_val;
      tag_size_val = t.tag_size_val;
      core = (pbam_core_32*)(read_buffer + 4);
      realized = true;
      validate();
    } else {
      reset();
    }    
  } else if(t.validate()) {
    read_buffer = t.read_buffer;
    block_size_val = t.block_size_val;
    tag_size_val = t.tag_size_val;
    core = (pbam_core_32*)(read_buffer + 4);
    realized = false;
    validate();
  } else {
    read_buffer = NULL;
    realized = false;
    core = NULL;
    block_size_val = 0;   tag_size_val = 0;
  }
}

// ************************** Copy assignment operator *************************
inline pbam1_t & pbam1_t::operator = (const pbam1_t &t)
{
  // Check for self assignment
  if(this != &t) {
    if(t.isReal()) {
      if(t.validate()) {
        read_buffer = (char*)malloc(t.block_size_val + 1);
        memcpy(read_buffer, t.read_buffer, t.block_size_val);
        block_size_val = t.block_size_val;
        tag_size_val = t.tag_size_val;
        core = (pbam_core_32*)(read_buffer + 4);
        realized = true;
        validate();
      } else {
        reset();
      }    
    } else if(t.validate()) {
      read_buffer = t.read_buffer;
      block_size_val = t.block_size_val;
      tag_size_val = t.tag_size_val;
      core = (pbam_core_32*)(read_buffer + 4);
      realized = false;
      validate();
    } else {
      read_buffer = NULL;
      realized = false;
      core = NULL;
      block_size_val = 0;   tag_size_val = 0;
    }
  }
  return *this;
}

#endif