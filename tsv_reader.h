#ifndef __TSV_READER_H
#define __TSV_READER_H

#include <cstdlib>
#include <cstring>
#include <climits>
#include <vector>

extern "C" {
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"
}
#include "Error.h"

// a class to read tab-limited (tabixable) file using htsFile.h and kstring.h
class tsv_reader {
public:
  std::string filename;  // file name to read
  htsFile* hp;           // handles to htsFile 
  tbx_t* tbx;            // tabix index, if exists
  hts_itr_t* itr;        // iterator, if exists
  kstring_t str;         // kstring_t object to store the string
  int32_t lstr;          // length of the string
  int32_t nfields;       // number of tokenized fields (from kstring_t)
  int32_t* fields;       // indices of the starting point to each field 
  int32_t nlines;        // total number of lines read
  int32_t delimiter;     // delimiter to tokenize
  
  bool open(const char* filename); // open a file, set up the file handle
  bool close();                    // close the file, returns false if fails
  int32_t read_line();             // read a line, returns the number of tokenzied fields (=nfields)
  const char* str_field_at(int32_t idx); // get a pointer to the string at index idx
  int32_t int_field_at(int32_t idx);     // get integer value at index idx
  double double_field_at(int32_t idx);   // get double value at index idx
  int32_t store_to_vector(std::vector<std::string>& v); // store the tokenized values into string vectors
  bool jump_to(const char* reg);   // jump to a specific region using tabix
  bool jump_to(const char* chr, int32_t beg, int32_t end = INT_MAX); // jump to a specific region using tabix

 tsv_reader() : hp(NULL), tbx(NULL), itr(NULL), lstr(0), nfields(0), fields(NULL), nlines(0), delimiter(0) {
    str.l = str.m = 0; str.s = NULL;
  }

 tsv_reader(const char* filename) : hp(NULL), tbx(NULL), itr(NULL), lstr(0), nfields(0), fields(NULL), nlines(0), delimiter(0) {
    str.l = str.m = 0; str.s = NULL;    
    if ( !open(filename) )
      error("[E:%s:%s %s] Cannot open file %s for reading", __FILE__, __LINE__, __FUNCTION__, filename);    
  }  

  ~tsv_reader() {
    if ( str.s ) free(str.s);
    if ( fields ) free(fields);
  }
};
#endif
