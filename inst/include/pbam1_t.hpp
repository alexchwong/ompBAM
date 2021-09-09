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
  private:
    // Variables
    char * read_buffer;
    bool realized = false;
    pbam_core_32 * core;
    uint32_t block_size_val; uint32_t tag_size_val;
    std::map< std::string, pbam_tag_index > tag_index;
    
    // Functions
    void reset();
    
    void seq_to_str(const uint8_t val, std::string & dest);
    char cigar_op_to_char(uint32_t cigar_op);
    void cigar_to_str(const uint32_t val, std::string & dest);
    
    void build_tag_index();
    
    char search_tag_type(const std::string tag);
    char search_tag_subtype(const std::string tag);
    uint32_t search_tag_pos(const std::string tag);
    uint32_t search_tag_length(const std::string tag);
  public:
    pbam1_t();
    ~pbam1_t();
    
    // Used by pbam_in::supplyRead()
    pbam1_t(char * src, const bool realize = false);       
    
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
    
    // Ask if pbam1_t is "real" (exist on separate buffer), or "virtual"
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
        
    // For long reads, the cigar is stored as a "CG" tag of type B,I
    // If "CG" tag exists, return its length; otherwise return n_cigar_op
    uint32_t cigar_size();
    
    /* 
      Buffer-based getters for variable-length data:
      - These functions return a pointer to the raw data
      - By avoiding assignment, some programs may be quicker
      - Be careful not to use these pointers to write to the read buffer
          otherwise you may corrupt the data.
    */
    char * read_name();                     // Direct char pointer
    uint32_t * cigar();                     // Direct uint32_t pointer
    uint8_t * seq();                        // Direct uint8_t pointer       
    char * qual();                          // Direct char pointer

    /*
      Reference-based getters for variable-length data:
      - These functions work by the user providing a reference to the variable
          to store the data.
      - These functions are likely more user-friendly as there is less formatting
          of the data involved.
      - They may consume an additional overhead due to the need to copy the
          data to another location
    */
    
    // Fills string reference with read name
    // - Returns l_read_name if success, or 0 if fail to validate
    uint8_t read_name(std::string & dest);   

    // Fills a string with the SAM-based cigar string
    // - Returns n_cigar_op if success
    // - Returns 0 if failed to validate
    int cigar(std::string & dest);

    // Returns the cigar operation / value (as char / uint32_t) 
    //   given the position of the cigar operation (0 <= pos < n_cigar_op)
    char cigar_op(const uint16_t pos);
    uint32_t cigar_val(const uint16_t pos);

    // Returns the cigar operations and values (as char / uint32_t) as a vector
    int cigar_ops_and_vals(std::vector<char> & ops, std::vector<uint32_t> & vals);

    // Fills string reference with read sequence
    // - Returns length of sequence if success, or zero if fail to validate
    int seq(std::string & dest);

    // Returns a vector of uint8_t of per-base quality scores
    // - returns l_seq if success or 0 if fail
    int qual(std::vector<uint8_t> & dest); 
    // From SAMv1.pdf:
    // Base qualities are stored as bytes in the range [0, 93], 
    //   without any +33 conversion to printable ASCII
    // When base qualities are omitted but the sequence is not, 
    //   qual is ﬁlled with 0xFF bytes (to length l_seq).
    // NB: 0xFF = 255
    
    /*
       Tag Getters:
       - Copies data contained within given tag to the dest buffer.
       - Returns -1 if the data type is not appropriate for given tag
       - Returns 0 if success
    */
    
    // Fills a vector with the tags available to each read
    int AvailTags(std::vector<std::string> & tags);

    // Returns a char of the size, type and subtype designates of the tag
    // - Tag_Size() returns 0 if the tag doesn't exist
    // - Tag_Type() returns '\0' if the tag doesn't exist
    // - Tag_Subtype() returns '\0' if the tag doesn't exist or is not of type 'B'
    char Tag_Type(const std::string tag);
    char Tag_Subtype(const std::string tag);
    uint32_t Tag_Size(const std::string tag);
    
    // Returns raw char pointer to the beginning of the info stored by the tag
    // - For advanced users
    char * p_tagVal(const std::string tag);

    // Returns values of fixed length
    // - For tags of type AcCsSiIf
    char tagVal_A(const std::string tag);
    int8_t tagVal_c(const std::string tag);
    uint8_t tagVal_C(const std::string tag);
    int16_t tagVal_s(const std::string tag);
    uint16_t tagVal_S(const std::string tag);
    int32_t tagVal_i(const std::string tag);
    uint32_t tagVal_I(const std::string tag);
    float tagVal_f(const std::string tag);

    // Returns a Z-tag by reference to a string
    int tagVal_Z(const std::string tag, std::string & dest);    // 'Z'

    // Returns a B-tag by reference to its respective type
    int tagVal_B(const std::string tag, std::vector<int8_t> & dest);   // 'B, c'
    int tagVal_B(const std::string tag, std::vector<uint8_t> & dest);    // 'B, C'
    int tagVal_B(const std::string tag, std::vector<int16_t> & dest);    // 'B, s'
    int tagVal_B(const std::string tag, std::vector<uint16_t> & dest);   // 'B, S'
    int tagVal_B(const std::string tag, std::vector<int32_t> & dest);    // 'B, i'
    int tagVal_B(const std::string tag, std::vector<uint32_t> & dest);   // 'B, I'
    int tagVal_B(const std::string tag, std::vector<float> & dest);      // 'B, f'
};

#include "pbam1_t_initializers.hpp"
#include "pbam1_t_getters.hpp"
#include "pbam1_t_tag_getters.hpp"
#include "pbam1_t_internals.hpp"


#endif