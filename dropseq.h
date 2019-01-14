#ifndef __DROPSEQ_H
#define __DROPSEQ_H

#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>
#include "Error.h"

#include "htslib/sam.h"

class DropCell;
class DropTranscript;

typedef std::map< int32_t, int32_t > ii_map_t;
typedef std::map< int32_t, int32_t >::iterator ii_map_it_t;
typedef std::pair< int32_t, int32_t > ii_pair_t;
typedef std::map< std::string, DropCell* > sc_map_t;
typedef std::map< std::string, DropCell* >::iterator sc_map_it_t;
typedef std::map< int32_t, DropTranscript* > it_map_t;
typedef std::map< int32_t, DropTranscript* >::iterator it_map_it_t;
typedef std::map< std::string, ii_pair_t > sr_map_t;
typedef std::map< std::string, ii_pair_t >::iterator sr_map_it_t;

// represent a dropseq library
class DropLibrary {
 public:
  sc_map_t mapCell;  // <string, DropCell*>
   
  ~DropLibrary();
  
  DropCell* getCell(std::string& barcode);
  ii_pair_t& addRead(std::string& barcode, int32_t tid, std::string& umi, int32_t nh);
  //ii_pair_t& addRead(bam_hdr_t* h, bam1_t* b, const char* NH, const char* NM);
};

// represent a dropseq cell
class DropCell {
 public:
  it_map_t mapTranscript; // <int, DropTranscript*>

  ~DropCell();

  DropTranscript* getTranscript(int32_t tid);
  ii_pair_t& addRead(int32_t tid, std::string& umi, int32_t nh);
};

class DropTranscript {
 public:
  sr_map_t mapRead;  // <string, <int,int> >  : UMI, <numDups, minNH>
  
  ii_pair_t& addRead( std::string& umi, int32_t nh);
};

#endif
