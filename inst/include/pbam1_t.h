#ifndef _pbam1_t
#define _pbam1_t

struct pbam_core_32{
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

class pbam1_t{
  private:
    char * read_buffer;
    bool realized;    // set to true if read_buffer is an allocated memory buffer (rather than just a pointer)
    
    pbam_core_32 * core;
    
    // Validation tests:
    uint32_t block_size;  // assigned on creation. check on validate()
    bool validate();
    void reset();
    
    void seq_to_str(uint8_t val, string & dest);
    
    uint32_t search_tag_pos(const string & tag);
  public:
    pbam1_t();
    pbam1_t(char * src, bool realize = false);
    ~pbam1_t();
    pbam1_t(const pbam1_t &t);
    pbam1_t & operator = (const pbam1_t &t);
    
    void realize();
    bool isReal() {return(realized);};
    
    // Getters:
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
    
    char * read_name();
    uint8_t read_name(std::string & str);
    
    uint32_t * cigar();
    uint16_t cigar(uint32_t * dest);
    
    uint8_t * seq();
    uint32_t seq(uint8_t * dest);
    int seq(string & dest);
    
    char * qual();
    uint32_t qual(char * dest);
    
    int tagVal(const string & tag, char * dest);
    int tagVal(const string & tag, int8_t * dest);
    int tagVal(const string & tag, uint8_t * dest);
    int tagVal(const string & tag, int16_t * dest);
    int tagVal(const string & tag, uint16_t * dest);
    int tagVal(const string & tag, int32_t * dest);
    int tagVal(const string & tag, uint32_t * dest);
    int tagVal(const string & tag, float * dest);
    int tagVal(const string & tag, string & dest);
    int tagVal(const string & tag, std::vector<int8_t> & dest);
    int tagVal(const string & tag, std::vector<uint8_t> dest);
    int tagVal(const string & tag, std::vector<int16_t> dest);
    int tagVal(const string & tag, std::vector<uint16_t> dest);
    int tagVal(const string & tag, std::vector<int32_t> dest);
    int tagVal(const string & tag, std::vector<uint32_t> dest);
    int tagVal(const string & tag, std::vector<float> dest);

};

// Constructor (empty):
pbam1_t::pbam1_t() {
  read_buffer = NULL;
  realized = false;
  core = NULL;
  block_size = 0;
}

void pbam1_t::reset() {
  if(read_buffer && realized) {
    free(read_buffer);
    read_buffer = NULL;
  }
  realized = false;
  core = NULL;
  block_size = 0;
}

// Constructor (takes char*, creates virtual read pointing to buffer):
pbam1_t::pbam1_t(char * src, bool realize) {
  uint32_t *temp_block_size = (uint32_t *)(src);
  block_size = *temp_block_size;
  if(realize) {
    read_buffer = (char*)malloc(block_size + 1);
    memcpy(read_buffer, src, block_size);
    core = (pbam_core_32 *)(read_buffer + 4);
    realized = true;
  } else {
    read_buffer = src;
    realized = false;
    core = (pbam_core_32 *)(read_buffer + 4);
  }
  validate();
}

// Destructor
pbam1_t::~pbam1_t() {
  if(read_buffer && realized) {
    free(read_buffer);
    read_buffer = NULL;
  }
  realized = false;
  core = NULL;
  block_size = 0;
}

void pbam1_t::realize() {
  if(validate() && !realized) {
    char *tmp = read_buffer;
    read_buffer = (char*)malloc(block_size + 1);
    memcpy(read_buffer, tmp, block_size);
    core = (pbam_core_32 *)(read_buffer + 4);
    realized = true;
  }
}

// copy constructor
pbam1_t(const pbam1_t &t) {
  if(t.isReal() && t.validate()) {
    read_buffer = (char*)malloc(t.block_size + 1);
    memcpy(read_buffer, t.read_buffer, t.block_size);
    block_size = t.block_size;
    core = (pbam_core_32*)(read_buffer + 4);
    realized = true;
  } else if(t.validate()) {
    read_buffer = t.read_buffer;
    block_size = t.block_size;
    core = (pbam_core_32*)(read_buffer + 4);
    realized = false;
  } else {
    read_buffer = NULL;
    realized = false;
    core = NULL;
    block_size = 0;
  }
}

// copy assignment operator
pbam1_t & pbam1_t::operator = (const pbam1_t &t)
{
  // Check for self assignment
  if(this != &t) {
    if(t.isReal() && t.validate()) {
      read_buffer = (char*)malloc(t.block_size + 1);
      memcpy(read_buffer, t.read_buffer, t.block_size);
      block_size = t.block_size;
      core = (pbam_core_32*)(read_buffer + 4);
      realized = true;
    } else if(t.validate()) {
      read_buffer = t.read_buffer;
      block_size = t.block_size;
      core = (pbam_core_32*)(read_buffer + 4);
      realized = false;
    } else {
      read_buffer = NULL;
      realized = false;
      core = NULL;
      block_size = 0;
    }
  }
  return *this;
}

// validate:
bool pbam1_t::validate() {
  if(block_size < 36) {
    reset(); return(false);
  }
  uint32_t *temp_block_size = (uint32_t *)(read_buffer);
  if(*temp_block_size != block_size) {
    reset(); return(false);
  }
  if(!core) {
    reset(); return(false);
  }
  if(
      (36 + core->l_read_name + 
        core->n_cigar_op * 4 + 
        core->l_seq + 
        ((core->l_seq + 1) / 2)
      ) > block_size
  ) {
    reset(); return(false);
  }
  return(true);
}

int32_t pbam1_t::refID() {
  if(validate()) {
    return(core->refID);
  }
  return(0);
}
int32_t pbam1_t::pos() {
  if(validate()) {
    return(core->pos);
  }
  return(0);
}
uint8_t pbam1_t::l_read_name() {
  if(validate()) {
    return(core->l_read_name);
  }
  return(0);
}
uint8_t pbam1_t::mapq(){
  if(validate()) {
    return(core->mapq);
  }
  return(0);
}
uint16_t pbam1_t::bin(){
  if(validate()) {
    return(core->bin);
  }
  return(0);
}
uint16_t pbam1_t::n_cigar_op(){
  if(validate()) {
    return(core->n_cigar_op);
  }
  return(0);
}
uint32_t pbam1_t::flag(){
  if(validate()) {
    return(core->flag);
  }
  return(0);
}
uint32_t pbam1_t::l_seq(){
  if(validate()) {
    return(core->l_seq);
  }
  return(0);
}
int32_t pbam1_t::next_refID(){
  if(validate()) {
    return(core->next_refID);
  }
  return(0);
}
int32_t pbam1_t::next_pos(){
  if(validate()) {
    return(core->next_pos);
  }
  return(0);
}
int32_t pbam1_t::tlen(){
  if(validate()) {
    return(core->tlen);
  }
  return(0);
}

char * pbam1_t::read_name() {
  if(validate()) {
    return((char*)(read_buffer + 36));
  }
  return(NULL);
}

uint8_t pbam1_t::read_name(std::string & str) {
  if(validate()) {
    char *tmp = (char*)(read_buffer + 36);
    str.assign(tmp);
    return(core->l_read_name);
  }
  return(-1);
}

// NB does not yet support long reads where cigar length > 65535
uint32_t * pbam1_t::cigar() {
  if(validate()) {
    return((uint32_t*)(read_buffer + 36 + core->l_read_name));
  }
  return(NULL);
}

uint16_t pbam1_t::cigar(uint32_t* dest) {
  if(validate()) {
    dest = (uint32_t*)malloc(sizeof(uint32_t) * (core->n_cigar_op + 1));
    memcpy(dest, read_buffer + 36 + core->l_read_name, sizeof(uint32_t) * (core->n_cigar_op));
    return(core->n_cigar_op);
  }
  return(0);
}

uint8_t * pbam1_t::seq() {
  if(validate()) {
    return((uint8_t*)(read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op));
  }
  return(NULL);
}

uint32_t pbam1_t::seq(uint8_t * dest) {
  if(validate()) {
    dest = (uint8_t*)malloc(sizeof(uint8_t) * ((core->l_seq + 1) / 2));
    memcpy(dest, read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op, 
      sizeof(uint8_t) * ((core->l_seq + 1) / 2));
    return(core->l_seq);
  }
  return(0);
}

void pbam1_t::seq_to_str(uint8_t val, string & dest) {
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

int pbam1_t::seq(string & dest) {
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

char * pbam1_t::qual() {
  if(validate()) {
    return((char*)(read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op + ((core->l_seq + 1) / 2)));
  }
  return(NULL);
}

uint32_t pbam1_t::seq(uint8_t * dest) {
  if(validate()) {
    dest = (uint8_t*)malloc(core->l_seq);
    memcpy(dest, read_buffer + 36 + core->l_read_name + sizeof(uint32_t) * core->n_cigar_op + ((core->l_seq + 1) / 2),
      core->l_seq);
    return(core->l_seq);
  }
  return(0);
}

// TAGS:
uint32_t search_tag_pos(const string & tag) {
  char null = '\0';
  uint32_t tag_pos = (36 + 
    core->l_read_name + 
    core->n_cigar_op * 4 + 
    core->l_seq + 
    ((core->l_seq + 1) / 2)
  );
  char *type, *subtype *byte; uint32_t *b_len;
  while(tag_pos < block_size) {
    if(strncmp(tag.c_str(), read_buffer + tag_pos, 2) == 0) return(tag_pos);
    type = read_buffer + tag_pos + 2;
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
        while(strncmp(null, byte, 1) != 0) {
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
          case 'S': case 'S':
            tag_pos += 8 + *b_len * 2; break;
          case 'i': case 'I': case 'f':
            tag_pos += 8 + *b_len * 4; break;
        }
    }
  }
  return(tag_pos);
}


int tagVal(const string & tag, char * dest) {
  if(validate()) {
    uint32_t tag_pos = search_tag_pos(tag);
  }
}

#endif