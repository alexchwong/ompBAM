
#ifndef _pbam1_t_tag_getters
#define _pbam1_t_tag_getters

// TAGS:

inline int pbam1_t::AvailTags(std::vector<std::string> & tags) {
  tags.clear();
  build_tag_index();
  for(std::map< std::string, pbam_tag_index >::iterator tagmap = tag_index.begin(); 
      tagmap != tag_index.end(); tagmap++) {
    tags.push_back(tagmap->first);
  }
  return(tag_index.size());
}

inline char pbam1_t::Tag_Type(const std::string tag) {
  return(search_tag_type(tag));
}

inline char pbam1_t::Tag_Subtype(const std::string tag) {
  return(search_tag_subtype(tag));
}

inline uint32_t pbam1_t::Tag_Size(const std::string tag) {
  return(search_tag_length(tag));
}

// Return a char pointer to the raw value contained in the tag
// - The user should typecast the pointer before obtaining the proper value
inline char * pbam1_t::p_tagVal(const std::string tag) {
  if(validate()) {
    char type = search_tag_type(tag);
    uint32_t tag_pos = search_tag_pos(tag);
    switch(type) {
      case 'c': case 'C': case 's': case 'S': case 'i': case 'I': case 'f': case 'Z':
        return(read_buffer + tag_pos + 3);
      case 'B':
        // User's responsibility to obtain length of this buffer
        return(read_buffer + tag_pos + 8);
      default:
        return(NULL);
    }
  }
  return(NULL);
}

inline char pbam1_t::tagVal_A(const std::string tag) {
  if(validate()) {
    if(search_tag_type(tag) == 'A') {
      uint32_t tag_pos = search_tag_pos(tag);
      char * tmp = read_buffer + tag_pos + 3;
      return(*tmp);
    }
  }
  return('\0');
}

inline int8_t pbam1_t::tagVal_c(const std::string tag) {
  if(validate()) {
    if(search_tag_type(tag) == 'c') {
      int8_t tag_pos = search_tag_pos(tag);
      char * tmp = read_buffer + tag_pos + 3;
      return(*(int8_t*)tmp);
    }
  }
  return(0);
}

inline uint8_t pbam1_t::tagVal_C(const std::string tag) {
  if(validate()) {
    if(search_tag_type(tag) == 'C') {
      uint32_t tag_pos = search_tag_pos(tag);
      char * tmp = read_buffer + tag_pos + 3;
      return(*(uint8_t*)tmp);
    }
  }
  return(0);
}

inline int16_t pbam1_t::tagVal_s(const std::string tag) {
  if(validate()) {
    if(search_tag_type(tag) == 's') {
      uint32_t tag_pos = search_tag_pos(tag);
      char * tmp = read_buffer + tag_pos + 3;
      return(*(int16_t*)tmp);
    }
  }
  return(0);
}


inline uint16_t pbam1_t::tagVal_S(const std::string tag) {
  if(validate()) {
    if(search_tag_type(tag) == 'S') {
      uint32_t tag_pos = search_tag_pos(tag);
      char * tmp = read_buffer + tag_pos + 3;
      return(*(uint16_t*)tmp);
    }
  }
  return(0);
}


inline int32_t pbam1_t::tagVal_i(const std::string tag) {
  if(validate()) {
    if(search_tag_type(tag) == 'i') {
      uint32_t tag_pos = search_tag_pos(tag);
      char * tmp = read_buffer + tag_pos + 3;
      return(*(int32_t*)tmp);
    }
  }
  return(0);
}


inline uint32_t pbam1_t::tagVal_I(const std::string tag) {
  if(validate()) {
    if(search_tag_type(tag) == 'I') {
      uint32_t tag_pos = search_tag_pos(tag);
      char * tmp = read_buffer + tag_pos + 3;
      return(*(uint32_t*)tmp);
    }
  }
  return(0);
}


inline float pbam1_t::tagVal_f(const std::string tag) {
  if(validate()) {
    if(search_tag_type(tag) == 'f') {
      uint32_t tag_pos = search_tag_pos(tag);
      char * tmp = read_buffer + tag_pos + 3;
      return(*(float*)tmp);
    }
  }
  return(0);
}


// For Z-type tags (contains a string of variable length)
inline int pbam1_t::tagVal_Z(const std::string tag, std::string & dest) {
  dest.clear();
  if(validate()) {
    char type = search_tag_type(tag);
    if(type == 'Z') {
      uint32_t tag_length = search_tag_length(tag);
      uint32_t tag_pos = search_tag_pos(tag);
      char * tmp = read_buffer + tag_pos + 3;
      dest.assign(std::string(tmp, tag_length));
      return(tag_length);
    }
  }
  return(-1);
}

// For B-type tags of type c
inline int pbam1_t::tagVal_B(const std::string tag, std::vector<int8_t> & dest) {
  dest.clear();
  if(validate()) {
    if(search_tag_type(tag) == 'B') {
      if(search_tag_subtype(tag) == 'c') {
        uint32_t tag_length = search_tag_length(tag);
        uint32_t tag_pos = search_tag_pos(tag);
        int8_t* tmp = (int8_t*)(read_buffer + tag_pos + 8);
        for(unsigned int i = 0; i < tag_length; i++) {
          dest.push_back(*tmp);
          tmp++;
        }
        return(tag_length);
      }
    }
  }
  return(-1);
}

inline int pbam1_t::tagVal_B(const std::string tag, std::vector<uint8_t> & dest) {
  dest.clear();
  if(validate()) {
    if(search_tag_type(tag) == 'B') {
      if(search_tag_subtype(tag) == 'C') {
        uint32_t tag_length = search_tag_length(tag);
        uint32_t tag_pos = search_tag_pos(tag);
        uint8_t* tmp = (uint8_t*)(read_buffer + tag_pos + 8);
        for(unsigned int i = 0; i < tag_length; i++) {
          dest.push_back(*tmp);
          tmp++;
        }
        return(tag_length);
      }
    }
  }
  return(-1);
}

inline int pbam1_t::tagVal_B(const std::string tag, std::vector<int16_t> & dest) {
  if(validate()) {
    dest.clear();
    if(search_tag_type(tag) == 'B') {
      if(search_tag_subtype(tag) == 's') {
        uint32_t tag_length = search_tag_length(tag);
        uint32_t tag_pos = search_tag_pos(tag);
        int16_t* tmp = (int16_t*)(read_buffer + tag_pos + 8);
        for(unsigned int i = 0; i < tag_length; i++) {
          dest.push_back(*tmp);
          tmp++;
        }
        return(tag_length);
      }
    }
  }
  return(-1);
}

inline int pbam1_t::tagVal_B(const std::string tag, std::vector<uint16_t> & dest) {
  if(validate()) {
    dest.clear();
    if(search_tag_type(tag) == 'B') {
      if(search_tag_subtype(tag) == 'S') {
        uint32_t tag_length = search_tag_length(tag);
        uint32_t tag_pos = search_tag_pos(tag);
        uint16_t* tmp = (uint16_t*)(read_buffer + tag_pos + 8);
        for(unsigned int i = 0; i < tag_length; i++) {
          dest.push_back(*tmp);
          tmp++;
        }
        return(tag_length);
      }
    }
  }
  return(-1);
}

inline int pbam1_t::tagVal_B(const std::string tag, std::vector<int32_t> & dest) {
  if(validate()) {
    dest.clear();
    if(search_tag_type(tag) == 'B') {
      if(search_tag_subtype(tag) == 'i') {
        uint32_t tag_length = search_tag_length(tag);
        uint32_t tag_pos = search_tag_pos(tag);
        int32_t* tmp = (int32_t*)(read_buffer + tag_pos + 8);
        for(unsigned int i = 0; i < tag_length; i++) {
          dest.push_back(*tmp);
          tmp++;
        }
        return(tag_length);
      }
    }
  }
  return(-1);
}

inline int pbam1_t::tagVal_B(const std::string tag, std::vector<uint32_t> & dest) {
  if(validate()) {
    dest.clear();
    if(search_tag_type(tag) == 'B') {
      if(search_tag_subtype(tag) == 'I') {
        uint32_t tag_length = search_tag_length(tag);
        uint32_t tag_pos = search_tag_pos(tag);
        uint32_t* tmp = (uint32_t*)(read_buffer + tag_pos + 8);
        for(unsigned int i = 0; i < tag_length; i++) {
          dest.push_back(*tmp);
          tmp++;
        }
        return(tag_length);
      }
    }
  }
  return(-1);
}

inline int pbam1_t::tagVal_B(const std::string tag, std::vector<float> & dest) {
  if(validate()) {
    dest.clear();
    if(search_tag_type(tag) == 'B') {
      if(search_tag_subtype(tag) == 'f') {
        uint32_t tag_length = search_tag_length(tag);
        uint32_t tag_pos = search_tag_pos(tag);
        float* tmp = (float*)(read_buffer + tag_pos + 8);
        for(unsigned int i = 0; i < tag_length; i++) {
          dest.push_back(*tmp);
          tmp++;
        }
        return(tag_length);
      }
    }
  }
  return(-1);
}

#endif