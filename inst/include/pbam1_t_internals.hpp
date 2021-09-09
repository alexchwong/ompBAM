#ifndef _pbam1_t_internals
#define _pbam1_t_internals

// ################################ INTERNAL FUNCTIONS #########################

// void seq_to_str(const uint8_t val, std::string & dest);
// void cigar_op_to_str(uint32_t cigar_op, std::string & dest);
// void cigar_to_str(const uint32_t val, std::string & dest);

// void build_tag_index();
// uint32_t search_tag_pos(const std::string tag, uint32_t & tag_length);
// uint32_t tag_check(const std::string tag, uint32_t & tag_length, char type);

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

// cigar_op ~ [0,8]
inline char pbam1_t::cigar_op_to_char(uint32_t cigar_op) {
  if(cigar_op <= 8) {
    switch(cigar_op) {
      case 0:
        return('M');
      case 1:
        return('I');
      case 2:
        return('D');
      case 3:
        return('N');
      case 4:
        return('S');
      case 5:
        return('H');
      case 6:
        return('P');
      case 7:
        return('=');
      case 8:
        return('X');
    }
  }
  return('\0');
}

inline void pbam1_t::cigar_to_str(uint32_t val, std::string & dest) {
  uint32_t len = val >> 4;
  uint32_t cigar_op = val & 15;
  if(cigar_op <= 8) {
    dest.append(std::to_string(len));
    dest.push_back(cigar_op_to_char(cigar_op));
  }
}

// The first time any tag is queried, run this to index every tag
inline void pbam1_t::build_tag_index() {  

  if(tag_index.size() == 0 && tag_size_val > 0) {
  // Ensures tag index is only built once and only when needed

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
      
      tag_index_entry.tag_pos = tag_pos;
      
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
    } // End of while loop to iterate over every tag
  }
}

inline char pbam1_t::search_tag_type(const std::string tag){
  if(tag_size_val == 0) return('\0');
  build_tag_index();
  auto it = tag_index.find(tag);
  if (it == tag_index.end()) return('\0');
  return(tag_index[tag].type);
}

inline char pbam1_t::search_tag_subtype(const std::string tag){
  if(tag_size_val == 0) return('\0');
  build_tag_index();
  auto it = tag_index.find(tag);
  if (it == tag_index.end()) return('\0');
  return(tag_index[tag].subtype);
}

inline uint32_t pbam1_t::search_tag_pos(const std::string tag) {
  if(tag_size_val == 0) return(0);
  build_tag_index();
  auto it = tag_index.find(tag);
  if (it == tag_index.end()) return(0);
  return(tag_index[tag].tag_pos);
}

inline uint32_t pbam1_t::search_tag_length(const std::string tag) {
  if(tag_size_val == 0) return(0);
  build_tag_index();
  auto it = tag_index.find(tag);
  if (it == tag_index.end()) return(0);
  return(tag_index[tag].tag_length);
}


#endif