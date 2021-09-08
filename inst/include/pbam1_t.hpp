#ifndef _pbam1_t
#define _pbam1_t

struct pbam_core_32{
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

// Tag type: acCsSiIfZB
// B Tag subtype: cCsSiIf
// Number of bytes from beginning of read_buffer
// Number of bytes from beginning of read_buffer
struct pbam_tag_index{
  char type;       
  char subtype;    
  uint32_t tag_pos;   
  uint32_t tag_length;   
};

class pbam1_t{
  public:
    pbam1_t();
    
    // Used by pbam_in::supplyRead()
    pbam1_t(char * src, const bool realize = false);    
    
    ~pbam1_t();
    
    // Copy constructor
    pbam1_t(const pbam1_t &t);

    // Copy assignment operator
    pbam1_t & operator = (const pbam1_t &t);  
    
    // Check if read is a valid data structure
    bool validate() const;

    // Call this to copy the read to a dedicated buffer. This will
    //   cause this read to be persistent and validate after
    //   further calls to pbam_in::fillReads()
    // When pbam_in::fillReads() is called then the information
    //   held by virtual reads is lost. This is because, by default,
    //   pbam_in::supplyRead() returns a pbam1_t that contains pointers
    //   directly to the data buffer (for quick memory access).
    int realize();
    bool isReal() const {return(realized);};
    
    // ******************************* Getters ********************************
    uint32_t block_size() {return(block_size_val);};
    
    // Core:
    int32_t refID();
    int32_t pos();
    uint8_t l_read_name();
    uint8_t mapq();
    uint16_t bin();
    uint16_t n_cigar_op();
    uint32_t flag();
    uint32_t l_seq();
    int32_t next_refID();
    int32_t next_pos();
    int32_t tlen();
    
    // Read name:
    char * read_name();                     // Direct char pointer
    uint8_t read_name(std::string & dest);   // Fills string reference with read name
    
    // Cigar:
    uint32_t * cigar();                     // Direct uint32_t pointer
    uint16_t cigar(uint32_t * dest);        // malloc/memcpy's buffer to dest
    int cigar(std::string & dest);
    uint32_t cigar_val(const uint16_t pos);
    int cigar_op(const uint16_t pos, std::string & dest);
    
    // Sequence:
    uint8_t * seq();                        // Direct uint8_t pointer       
    uint32_t seq(uint8_t * dest);           // malloc/memcpy's buffer to dest
    int seq(std::string & dest);            // Fills string reference with sequence
                                            //   returns l_seq if success
    
    char * qual();                          // Direct char pointer
    uint32_t qual(uint8_t * dest);          // memcpy's converted quality scores; 
                                            //   returns l_seq if success
    // From SAMv1.pdf:
    // Base qualities are stored as bytes in the range [0, 93], 
    //   without any +33 conversion to printable ASCII
    // When base qualities are omitted but the sequence is not, 
    //   qual is ﬁlled with 0xFF bytes (to length l seq).
    
    // Tag Getters:
    // - Copies data contained within given tag to the dest buffer.
    // - Returns -1 if the data type is not appropriate for given tag
    // - Returns 0 if success
    int AvailTags(std::vector<std::string> & tags);
    
    int tagVal(const std::string tag, char * dest);      // 'a'
    int tagVal(const std::string tag, int8_t * dest);    // 'c'
    int tagVal(const std::string tag, uint8_t * dest);   // 'C'
    int tagVal(const std::string tag, int16_t * dest);   // 's'
    int tagVal(const std::string tag, uint16_t * dest);  // 'S'
    int tagVal(const std::string tag, int32_t * dest);   // 'i'
    int tagVal(const std::string tag, uint32_t * dest);  // 'I'
    int tagVal(const std::string tag, float * dest);     // 'f'
    int tagVal(const std::string tag, std::string & dest);    // 'Z'
    int tagVal(const std::string tag, std::vector<int8_t> & dest);   // 'B, c'
    int tagVal(const std::string tag, std::vector<uint8_t> & dest);    // 'B, C'
    int tagVal(const std::string tag, std::vector<int16_t> & dest);    // 'B, s'
    int tagVal(const std::string tag, std::vector<uint16_t> & dest);   // 'B, S'
    int tagVal(const std::string tag, std::vector<int32_t> & dest);    // 'B, i'
    int tagVal(const std::string tag, std::vector<uint32_t> & dest);   // 'B, I'
    int tagVal(const std::string tag, std::vector<float> & dest);      // 'B, f'

  private:
    char * read_buffer;
    bool realized;    // set to true if read_buffer is an 
                      // allocated memory buffer (rather than just a pointer)
    pbam_core_32 * core;
    uint32_t block_size_val;  // assigned on creation. check on validate()
    uint32_t tag_size_val;    // assigned on creation. check on validate()  
    
    std::map< std::string, pbam_tag_index > tag_index;
    
    void reset();
    
    void seq_to_str(const uint8_t val, std::string & dest);
    void cigar_op_to_str(uint32_t cigar_op, std::string & dest);
    void cigar_to_str(const uint32_t val, std::string & dest);
    
    void build_tag_index();
    uint32_t search_tag_pos(const std::string tag, uint32_t & tag_length);
    uint32_t tag_check(const std::string tag, uint32_t & tag_length, char type);
};

// ************************** Constructor (empty) ******************************
inline pbam1_t::pbam1_t() {
  read_buffer = NULL;
  realized = false;
  core = NULL;
  block_size_val = 0;   tag_size_val = 0;
}

// ******** Constructor given pointer to buffer; option to realize read ********
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

// ************************ Copies read to real buffer ************************
inline int pbam1_t::realize() {
  if(validate() && !realized) {
    char *tmp = read_buffer;
    read_buffer = (char*)malloc(block_size_val + 1);
    memcpy(read_buffer, tmp, block_size_val);
    core = (pbam_core_32 *)(read_buffer + 4);
    realized = true;
  }
  if(validate()) {
    return(0);
  }
  return(-1);
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

// ********************************* Validity *********************************

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
    Rcpp::Rcout << "Invalid read: " << read_name_s << "\n";
    return(false);
  }
  
  return(true);
}

// ********************************* Getters **********************************

inline int32_t pbam1_t::refID() {
  if(validate()) {
    return(core->refID);
  }
  return(0);
}
inline int32_t pbam1_t::pos() {
  if(validate()) {
    return(core->pos);
  }
  return(0);
}
inline uint8_t pbam1_t::l_read_name() {
  if(validate()) {
    return(core->l_read_name);
  }
  return(0);
}
inline uint8_t pbam1_t::mapq(){
  if(validate()) {
    return(core->mapq);
  }
  return(0);
}
inline uint16_t pbam1_t::bin(){
  if(validate()) {
    return(core->bin);
  }
  return(0);
}
inline uint16_t pbam1_t::n_cigar_op(){
  if(validate()) {
    return(core->n_cigar_op);
  }
  return(0);
}
inline uint32_t pbam1_t::flag(){
  if(validate()) {
    return(core->flag);
  }
  return(0);
}
inline uint32_t pbam1_t::l_seq(){
  if(validate()) {
    return(core->l_seq);
  }
  return(0);
}
inline int32_t pbam1_t::next_refID(){
  if(validate()) {
    return(core->next_refID);
  }
  return(0);
}
inline int32_t pbam1_t::next_pos(){
  if(validate()) {
    return(core->next_pos);
  }
  return(0);
}
inline int32_t pbam1_t::tlen(){
  if(validate()) {
    return(core->tlen);
  }
  return(0);
}

inline char * pbam1_t::read_name() {
  if(validate()) {
    return((char*)(read_buffer + 36));
  }
  return(NULL);
}

inline uint8_t pbam1_t::read_name(std::string & dest) {
  if(validate()) {
    char *tmp = (char*)(read_buffer + 36);
    dest.assign(tmp);
    return(core->l_read_name);
  }
  return(-1);
}

// NB does not yet support long reads where cigar length > 65535
inline uint32_t * pbam1_t::cigar() {
  if(validate()) {
    return((uint32_t*)(read_buffer + 36 + core->l_read_name));
  }
  return(NULL);
}

inline uint16_t pbam1_t::cigar(uint32_t* dest) {
  if(validate()) {
    dest = (uint32_t*)malloc(sizeof(uint32_t) * (core->n_cigar_op + 1));
    memcpy(dest, read_buffer + 36 + core->l_read_name, sizeof(uint32_t) * (core->n_cigar_op));
    return(core->n_cigar_op);
  }
  return(0);
}

inline uint8_t * pbam1_t::seq() {
  if(validate()) {
    return((uint8_t*)(read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op));
  }
  return(NULL);
}

inline uint32_t pbam1_t::seq(uint8_t * dest) {
  if(validate()) {
    dest = (uint8_t*)malloc(sizeof(uint8_t) * ((core->l_seq + 1) / 2));
    memcpy(dest, read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op, 
      sizeof(uint8_t) * ((core->l_seq + 1) / 2));
    return(core->l_seq);
  }
  return(0);
}

inline void pbam1_t::cigar_op_to_str(uint32_t cigar_op, std::string & dest) {
  if(cigar_op <= 8) {
    switch(cigar_op) {
      case 0:
        dest.append("M"); break;
      case 1:
        dest.append("I"); break;
      case 2:
        dest.append("D"); break;
      case 3:
        dest.append("N"); break;
      case 4:
        dest.append("S"); break;
      case 5:
        dest.append("H"); break;
      case 6:
        dest.append("P"); break;
      case 7:
        dest.append("="); break;
      case 8:
        dest.append("X"); break; 
    }
  }
}

inline void pbam1_t::cigar_to_str(uint32_t val, std::string & dest) {
  uint32_t len = val >> 4;
  uint32_t cigar_op = val & 15;
  if(cigar_op <= 8) {
    dest.append(std::to_string(len));
    std::string cigar_op_str;
    cigar_op_to_str(cigar_op, cigar_op_str);
    dest.append(cigar_op_str);
  }
}

inline void pbam1_t::seq_to_str(uint8_t val, std::string & dest) {
  switch(val) {
    case 1:
      dest.append("A"); break;
    case 2:
      dest.append("C"); break;
    case 4:
      dest.append("G"); break;
    case 8:
      dest.append("T"); break;
    case 15:
      dest.append("N"); break;
    case 3:
      dest.append("M"); break;
    case 5:
      dest.append("R"); break;
    case 6:
      dest.append("S"); break;
    case 7:
      dest.append("V"); break;
    case 9:
      dest.append("W"); break;
    case 10:
      dest.append("Y"); break;
    case 11:
      dest.append("H"); break;
    case 12:
      dest.append("K"); break;
    case 13:
      dest.append("D"); break;
    case 14:
      dest.append("B"); break;
    case 0:
      dest.append("="); break;        
  }
}

inline int pbam1_t::seq(std::string & dest) {
  if(validate()) {
    dest.clear();
    uint8_t *tmp_buffer = (uint8_t *)(read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op);
    uint8_t *tmp, val;
    for(unsigned int i = 0; i < core->l_seq; i++) {
      tmp = tmp_buffer + (i/2);
      if(i % 2 == 0) {
        // odd
        val = *tmp >> 4;
      } else {
        // even
        val = *tmp % 16;
      }
      seq_to_str(val, dest);
    }
    return(core->l_seq);
  }
  return(0);
}

inline int pbam1_t::cigar(std::string & dest) {
  if(validate()) {
    dest.clear();
    char *tmp_buffer = (read_buffer + 36 + core->l_read_name);
    uint32_t *tmp, val;
    for(unsigned int i = 0; i < core->n_cigar_op; i++) {
      tmp = (uint32_t *)(tmp_buffer + (i * sizeof(uint32_t)));
      val = *tmp;
      cigar_to_str(val, dest);
    }
    return(core->n_cigar_op);
  }
  return(0);
}

inline uint32_t pbam1_t::cigar_val(const uint16_t pos) {
  if(pos < core->n_cigar_op) {
    char *tmp_buffer = (read_buffer + 36 + core->l_read_name);
    uint32_t *tmp = (uint32_t *)(tmp_buffer + (pos * sizeof(uint32_t)));
    return(*tmp >> 4);
  } else {
    return(0);
  }
}

inline int pbam1_t::cigar_op(const uint16_t pos, std::string & dest) {
  if(pos < core->n_cigar_op) {
    char *tmp_buffer = (read_buffer + 36 + core->l_read_name);
    uint32_t *tmp = (uint32_t *)(tmp_buffer + (pos * sizeof(uint32_t)));
    dest.clear();
    cigar_op_to_str(*tmp, dest);
    return(0);
  } else {
    return(-1);
  }
}

inline char * pbam1_t::qual() {
  if(validate()) {
    return((char*)(read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op + ((core->l_seq + 1) / 2)));
  }
  return(NULL);
}

inline uint32_t pbam1_t::qual(uint8_t * dest) {
  if(validate()) {
    dest = (uint8_t*)malloc(core->l_seq);
    memcpy(dest, read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op + ((core->l_seq + 1) / 2),
      core->l_seq);
    return(core->l_seq);
  }
  return(0);
}

// TAGS:

inline void pbam1_t::build_tag_index() {
  // Ensures tag index is only built once and only when needed
  if(tag_index.size() == 0 && tag_size_val > 0) {
    char null = '\0';
    uint32_t tag_pos = (36 + 
      core->l_read_name + 
      core->n_cigar_op * 4 + 
      ((core->l_seq + 1) / 2) +
      core->l_seq
    );
    
    pbam_tag_index tag_index_entry;
    char *type, *subtype, *byte; 
    uint32_t *b_len, tag_length;
    
    while(tag_pos < block_size_val + 4) {
      std::string tag_name(read_buffer + tag_pos, 2);
      
      type = read_buffer + tag_pos + 2;
      tag_index_entry.type = *type;
      
      if(*type == 'B') {
        subtype = read_buffer + tag_pos + 3;
        tag_index_entry.subtype = *subtype;
      }

      // Calculate tag length
      tag_length = 0;
      switch(*type) {
        case 'a': case 'c': case 'C': case 's': case 'S': case 'i': case 'I': case 'f':
          tag_length = 1; break;
        case 'Z':
          tag_length = 1;
          byte = read_buffer + tag_pos + 2 + tag_length;
          while(strncmp(&null, byte, 1) != 0) {
            tag_length++; 
            byte = read_buffer + tag_pos + 2 + tag_length;
          }
          break;
        case 'B':
          b_len = (uint32_t*)(read_buffer + tag_pos + 4);
          tag_length = *b_len; break;
      }

      tag_index_entry.tag_length = tag_length;
      
      tag_index.insert({tag_name, tag_index_entry});

      // Increment tag_pos for next entry:
      switch(*type) {
        case 'a': case 'c': case 'C':
          tag_pos += 4; break;
        case 's': case 'S':
          tag_pos += 5; break;
        case 'i': case 'I': case 'f':
          tag_pos += 7; break;
        case 'Z':
          tag_pos += 3;
          byte = read_buffer + tag_pos;
          while(strncmp(&null, byte, 1) != 0) {
            tag_pos++; 
            byte = read_buffer + tag_pos;
          }
          tag_pos++; break;
        case 'B':
          subtype = read_buffer + tag_pos + 3;
          b_len = (uint32_t*)(read_buffer + tag_pos + 4);
          switch(*subtype) {
            case 'c': case 'C':
              tag_pos += 8 + *b_len; break;
            case 's': case 'S':
              tag_pos += 8 + *b_len * 2; break;
            case 'i': case 'I': case 'f':
              tag_pos += 8 + *b_len * 4; break;
          }
          break;
        default:
          tag_pos = block_size_val + 4;
          Rcpp::Rcout << "Tag error - type not defined\n";
      }
    }
  }
}

inline uint32_t pbam1_t::search_tag_pos(const std::string tag, uint32_t & tag_length) {
  if(tag_size_val == 0) return(block_size_val + 4);
  build_tag_index();
  pbam_tag_index tag_index_entry = tag_index[tag];
  tag_length = tag_index_entry.tag_length;
  return(tag_index_entry.tag_pos);
}

inline uint32_t pbam1_t::tag_check(const std::string tag, uint32_t & tag_length, char type) {
  uint32_t tag_pos = search_tag_pos(tag, tag_length);
  if(tag_pos >= block_size_val + 4) return(0);
  char *tag_type = read_buffer + tag_pos + 2;
  char B = 'B';
  if(strncmp(tag_type, &B, 1) == 0) tag_type = read_buffer + tag_pos + 3;
  if(strncmp(tag_type, &type, 1) != 0) return(0);
  return(tag_pos);
}

inline int pbam1_t::AvailTags(std::vector<std::string> & tags) {
  tags.clear();
  build_tag_index();
  for(auto tagmap = tag_index.begin(); tagmap != tag_index.end(); tagmap++) {
    tags.push_back(tagmap->first);
  }
  return(tag_index.size());
}

inline int pbam1_t::tagVal(const std::string tag, char * dest) {
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'a');
    if(tag_pos == 0) return(-1);
    memcpy(dest, read_buffer + tag_pos + 3, 1);
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, int8_t * dest) {
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'c');
    if(tag_pos == 0) return(-1);
    memcpy(dest, read_buffer + tag_pos + 3, 1);
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, uint8_t * dest) {
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'C');
    if(tag_pos == 0) return(-1);
    memcpy(dest, read_buffer + tag_pos + 3, 1);
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, int16_t * dest) {
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 's');
    if(tag_pos == 0) return(-1);
    memcpy(dest, read_buffer + tag_pos + 3, 1);
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, uint16_t * dest) {
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'S');
    if(tag_pos == 0) return(-1);
    memcpy(dest, read_buffer + tag_pos + 3, 1);
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, int32_t * dest) {
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'i');
    if(tag_pos == 0) return(-1);
    memcpy(dest, read_buffer + tag_pos + 3, 1);
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, uint32_t * dest) {
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'I');
    if(tag_pos == 0) return(-1);
    memcpy(dest, read_buffer + tag_pos + 3, 1);
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, float * dest) {
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'f');
    if(tag_pos == 0) return(-1);
    memcpy(dest, read_buffer + tag_pos + 3, 1);
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, std::string & dest) {
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'Z');
    if(tag_pos == 0) return(-1);
    char * dest_c = (char*)malloc(tag_length + 1);
    memcpy(dest_c, read_buffer + tag_pos + 3, tag_length);
    dest.assign(std::string(dest_c, tag_length));
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, std::vector<int8_t> & dest){
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'c');
    if(tag_pos == 0) return(-1);
    
    char B = 'B';
    int8_t * dest_single = (int8_t*)malloc(sizeof(int8_t));
    dest.clear();
    
    if(strncmp(read_buffer + tag_pos + 2, &B, 1) == 0) {
      for(uint32_t i = 0; i < tag_length; i++) {
        memcpy(dest_single, read_buffer + tag_pos + 8 + i, sizeof(int8_t));
        dest.push_back(*dest_single);
      }
    } else {
      memcpy(dest_single, read_buffer + tag_pos + 3, sizeof(int8_t));
      dest.push_back(*dest_single);
    }
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, std::vector<uint8_t> & dest){
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'C');
    if(tag_pos == 0) return(-1);
    
    char B = 'B';
    uint8_t * dest_single = (uint8_t*)malloc(sizeof(uint8_t));
    dest.clear();
    
    if(strncmp(read_buffer + tag_pos + 2, &B, 1) == 0) {
      for(uint32_t i = 0; i < tag_length; i++) {
        memcpy(dest_single, read_buffer + tag_pos + 8 + i, sizeof(uint8_t));
        dest.push_back(*dest_single);
      }
    } else {
      memcpy(dest_single, read_buffer + tag_pos + 3, sizeof(uint8_t));
      dest.push_back(*dest_single);
    }
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, std::vector<int16_t> & dest){
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 's');
    if(tag_pos == 0) return(-1);
    
    char B = 'B';
    int16_t * dest_single = (int16_t*)malloc(sizeof(int16_t));
    dest.clear();
    
    if(strncmp(read_buffer + tag_pos + 2, &B, 1) == 0) {
      for(uint32_t i = 0; i < tag_length; i++) {
        memcpy(dest_single, read_buffer + tag_pos + 8 + i, sizeof(int16_t));
        dest.push_back(*dest_single);
      }
    } else {
      memcpy(dest_single, read_buffer + tag_pos + 3, sizeof(int16_t));
      dest.push_back(*dest_single);
    }
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, std::vector<uint16_t> & dest){
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'S');
    if(tag_pos == 0) return(-1);
    
    char B = 'B';
    uint16_t * dest_single = (uint16_t*)malloc(sizeof(uint16_t));
    dest.clear();
    
    if(strncmp(read_buffer + tag_pos + 2, &B, 1) == 0) {
      for(uint32_t i = 0; i < tag_length; i++) {
        memcpy(dest_single, read_buffer + tag_pos + 8 + i, sizeof(uint16_t));
        dest.push_back(*dest_single);
      }
    } else {
      memcpy(dest_single, read_buffer + tag_pos + 3, sizeof(uint16_t));
      dest.push_back(*dest_single);
    }
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, std::vector<int32_t> & dest){
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'i');
    if(tag_pos == 0) return(-1);
    
    char B = 'B';
    int32_t * dest_single = (int32_t*)malloc(sizeof(int32_t));
    dest.clear();
    
    if(strncmp(read_buffer + tag_pos + 2, &B, 1) == 0) {
      for(uint32_t i = 0; i < tag_length; i++) {
        memcpy(dest_single, read_buffer + tag_pos + 8 + i, sizeof(int32_t));
        dest.push_back(*dest_single);
      }
    } else {
      memcpy(dest_single, read_buffer + tag_pos + 3, sizeof(int32_t));
      dest.push_back(*dest_single);
    }
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, std::vector<uint32_t> & dest){
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'I');
    if(tag_pos == 0) return(-1);
    
    char B = 'B';
    uint32_t * dest_single = (uint32_t*)malloc(sizeof(uint32_t));
    dest.clear();
    
    if(strncmp(read_buffer + tag_pos + 2, &B, 1) == 0) {
      for(uint32_t i = 0; i < tag_length; i++) {
        memcpy(dest_single, read_buffer + tag_pos + 8 + i, sizeof(uint32_t));
        dest.push_back(*dest_single);
      }
    } else {
      memcpy(dest_single, read_buffer + tag_pos + 3, sizeof(uint32_t));
      dest.push_back(*dest_single);
    }
    return(0);
  }
  return(-1);
}

inline int pbam1_t::tagVal(const std::string tag, std::vector<float> & dest){
  if(validate()) {
    uint32_t tag_length;
    uint32_t tag_pos = tag_check(tag, tag_length, 'f');
    if(tag_pos == 0) return(-1);
    
    char B = 'B';
    float * dest_single = (float*)malloc(sizeof(float));
    dest.clear();
    
    if(strncmp(read_buffer + tag_pos + 2, &B, 1) == 0) {
      for(uint32_t i = 0; i < tag_length; i++) {
        memcpy(dest_single, read_buffer + tag_pos + 8 + i, sizeof(float));
        dest.push_back(*dest_single);
      }
    } else {
      memcpy(dest_single, read_buffer + tag_pos + 3, sizeof(float));
      dest.push_back(*dest_single);
    }
    return(0);
  }
  return(-1);
}

#endif