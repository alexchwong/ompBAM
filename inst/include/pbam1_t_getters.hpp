
#ifndef _pbam1_t_getters
#define _pbam1_t_getters

// ***************************** Core Getters **********************************

// These are easy to construct since they return a fixed-size return value

inline int32_t pbam1_t::refID() {
  if(validate()) return(core->refID);
  return(0);
}
inline int32_t pbam1_t::pos() {
  if(validate()) return(core->pos);
  return(0);
}
inline uint8_t pbam1_t::l_read_name() {
  if(validate()) return(core->l_read_name);
  return(0);
}
inline uint8_t pbam1_t::mapq(){
  if(validate()) return(core->mapq);
  return(0);
}
inline uint16_t pbam1_t::bin(){
  if(validate()) return(core->bin);
  return(0);
}
inline uint16_t pbam1_t::n_cigar_op(){
  if(validate()) return(core->n_cigar_op);
  return(0);
}
inline uint32_t pbam1_t::flag(){
  if(validate()) return(core->flag);
  return(0);
}
inline uint32_t pbam1_t::l_seq(){
  if(validate()) return(core->l_seq);
  return(0);
}
inline int32_t pbam1_t::next_refID(){
  if(validate()) return(core->next_refID);
  return(0);
}
inline int32_t pbam1_t::next_pos(){
  if(validate()) return(core->next_pos);
  return(0);
}
inline int32_t pbam1_t::tlen(){
  if(validate()) return(core->tlen);
  return(0);
}

// This is to accomodate for long reads which may contain >65535 operations
// For long reads, the cigar is stored as a "CG" tag of type B,I
// If "CG" tag exists, return its length; otherwise return n_cigar_op
inline uint32_t pbam1_t::cigar_size() {
  if(!validate()) return(0);
  if(core->n_cigar_op == 2) {
    uint32_t* c1 = (uint32_t*)(read_buffer + 36 + core->l_read_name)
    uint32_t* c2 = (uint32_t*)(read_buffer + 40 + core->l_read_name)
    if(
      cigar_op_to_char(*c1) == 'S' &&
      cigar_op_to_char(*c2) == 'N' &&
      *c1 >> 4 == core->l_seq
    ) {
      uint32_t size = search_tag_length("CG");
      if(size > 65535) return(size);
      return(0);
    }
  }
  return((uint32_t)core->n_cigar_op);
}

// ************************** Buffer-based  Getters ****************************
/* 
  Buffer-based getters for variable-length data:
  - These functions return a pointer to the raw data
  - By avoiding assignment, some programs may be quicker
  - Be careful not to use these pointers to write to the read buffer
      otherwise you may corrupt the data.
*/

inline char * pbam1_t::read_name() {
  if(validate()) return((char*)(read_buffer + 36));
  return(NULL);
}

inline uint32_t * pbam1_t::cigar() {
  if(validate()) {
    if(cigar_size() > 65535) return((uint32_t*)p_tagVal("CG"));
    return((uint32_t*)(read_buffer + 36 + core->l_read_name));
  }
  return(NULL);
}

inline uint8_t * pbam1_t::seq() {
  if(validate()) return((uint8_t*)(
    read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op));
  return(NULL);
}

inline char * pbam1_t::qual() {
  if(validate()) return((char*)(
    read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op + 
      ((core->l_seq + 1) / 2)));
  return(NULL);
}

/*
  Reference-based getters for variable-length data:
  - These functions work by the user providing a reference to the variable
      to store the data.
  - These functions are likely more user-friendly as there is less formatting
      of the data involved.
  - They may consume an additional overhead due to the need to copy the
      data to another location
*/

inline uint8_t pbam1_t::read_name(std::string & dest) {
  dest.clear();
  if(!validate()) return(0);
  char *tmp = (char*)(read_buffer + 36);
  dest.assign(tmp);
  return(core->l_read_name);
}

inline int pbam1_t::cigar(std::string & dest) {
  dest.clear();
  if(!validate()) return(0);
  uint32_t *tmp = cigar();
  uint32_t size = cigar_size();
  for(unsigned int i = 0; i < size; i++) {
    cigar_to_str(*tmp, dest);
    tmp++;
  }
  return(size);
}

inline char pbam1_t::cigar_op(const uint16_t pos) {
  if(!validate()) return('\0');
  if(pos < cigar_size()) {
    uint32_t *tmp = cigar() + pos;
    return(cigar_op_to_char(*tmp));
  }
  return('\0');
}

inline uint32_t pbam1_t::cigar_val(const uint16_t pos) {
  if(!validate()) return(0);
  if(pos < core->n_cigar_op) {
    uint32_t *tmp = cigar() + pos;
    return(*tmp >> 4);
  }
  return(0);
}

inline int pbam1_t::cigar_ops_and_vals(
  std::vector<char> & ops, std::vector<uint32_t> & vals
) {
  ops.clear(); vals.clear();
  if(!validate()) return(0);
  uint32_t size = cigar_size();
  if(size == 0) return(0);

  uint32_t *tmp = cigar();
  for(uint32_t i = 0; i < size; i++) {
    ops.push_back(cigar_op_to_char(*tmp));
    vals.push_back(*tmp >> 4);
    tmp++;
  }
  return(size);
}

inline int pbam1_t::seq(std::string & dest) {
  dest.clear();
  if(!validate())  return(0);
  
  uint8_t *tmp = (uint8_t *)(
    read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op);
  uint8_t val;
  for(unsigned int i = 0; i < core->l_seq; i++) {
    if(i % 2 == 0) {
      val = *tmp >> 4;
    } else {
      val = *tmp % 16;
      tmp++;
    }
    seq_to_str(val, dest);
  }
  return(core->l_seq);
}

inline int pbam1_t::qual(std::vector<uint8_t> & dest) {
  dest.clear();
  if(!validate()) return(0);
  uint8_t *tmp = (uint8_t*)(
    read_buffer + 36 + core->l_read_name + 
      (sizeof(uint32_t) * core->n_cigar_op) + 
      ((core->l_seq + 1) / 2)
  );
  for(unsigned int i = 0; i < core->l_seq; i++) {
    dest.push_back(*tmp);
    tmp++;
  }
  return(core->l_seq);
}

#endif