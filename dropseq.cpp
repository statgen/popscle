#include <cstring>

#include "dropseq.h"


// delete dropLibrary objects, first by deleting all attached cells.
DropLibrary::~DropLibrary() {
  for(sc_map_it_t it = mapCell.begin(); it != mapCell.end(); ++it) {
    delete it->second;
  }
}

// get cell info based on barcode, and create one if one does not exist
DropCell* DropLibrary::getCell(std::string& barcode) {
  sc_map_it_t it = mapCell.find(barcode);
  
  if ( it == mapCell.end() ) {
    return ( mapCell[barcode] = new DropCell() );
  }
  else {
    return it->second;
  }
}

// add a read into the library
ii_pair_t& DropLibrary::addRead(std::string& barcode, int32_t tid, std::string& umi, int32_t nh) {
  if ( rand() % 10000000 == 0 ) {
    notice("DropLibrary::addRead() - printing randomly every 10M reads mapCell.size() = %d",mapCell.size());
  }
  return getCell(barcode)->addRead(tid, umi, nh);
}

/*
ii_pair_t& DropLibrary::addRead(bam_hdr_t* h, bam1_t* b, const char* NH, const char* NM) {
  uint8_t* s = bam_aux_get(b, NH);
  if ( !s ) {
    error("Cannot find %c%c tag in record\n", NH[0], NH[1]);
  }
  int32_t vNH = bam_aux2i(s);

  s = bam_aux_get(b,NM);
  if ( !s ) {
    error("Cannot find %c%c tag in record\n", NM[0], NM[1]);
  }
  int vNM = bam_aux2i(s);

  // extract barcode
  char *prn = bam_get_qname(b);
  char *pbc = NULL;
  char *pumi = prn;
  char *ptmp = NULL;
  while( ( ptmp = strchr(pumi, ':') ) != NULL ) {
    pbc = pumi+1;
    pumi = ptmp+1;
  }

  return addRead(std::string(pbc,pumi-pbc-1), b->core.tid, std::string(pumi), vNH);  
}
*/

// delete dropCell objects, first by deleting all attached transcripts
DropCell::~DropCell() {
  for(it_map_it_t it = mapTranscript.begin(); it != mapTranscript.end(); ++it) {
    delete it->second;
  }
}

// get transcript info based on transcript ID, and create a new one if the transcript does not exist
DropTranscript* DropCell::getTranscript(int32_t tid) {
  it_map_it_t it = mapTranscript.find(tid);

  if ( it == mapTranscript.end() ) {
    return ( mapTranscript[tid] = new DropTranscript() );
  }
  else {
    return it->second;
  }
}

// add a read into the cell
ii_pair_t& DropCell::addRead(int32_t tid, std::string& umi, int32_t nh) {
  return getTranscript(tid)->addRead(umi, nh);
}

// add a read into the transcript
ii_pair_t& DropTranscript::addRead(std::string& umi, int32_t nh) {
  sr_map_it_t it = mapRead.find(umi);
  if ( it == mapRead.end() ) {
    return (mapRead[umi] = std::make_pair(1, nh));
  }
  else {
    ii_pair_t& ii = it->second;
    ++(ii.first);            // count 
    if ( ii.second > nh )    // retain mininum counts
      ii.second = nh;
    return ii;
  }
}


