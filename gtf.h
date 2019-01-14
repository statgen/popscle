#ifndef __APIGENOME_GTF_H
#define __APIGENOME_GTF_H

/////////////////////
// GTF format
// 1. seqname
// 2. source
// 3. feature - gene, transcript, exon, cds
// 4. start
// 5. end
// 6. score
// 7. strand
// 8. frame
// 9. attribute
////////////////////


/////////////////////////////////////////////////////////////////////////
/// GTF class should be able to (minimally)
/// 1. Read a GTF file and organize a structure
/// 2. Efficiently find genes and regions overlapping with certain position
///////////////////////////////////////////////////////////////////////////

#include <map>
#include <set>
#include <vector>
#include <cstring>
#include <climits>
#include "Error.h"
#include "gtf_interval_tree.h"

// A single genomic int32_t interval
class posLocus {
 public:
  int32_t beg1; // includes 1-based, excludes 0-based
  int32_t end0; // excludes 1-based, includes 0-based

  // constructor
  // b : beg1 value
  // e : end0 value
  posLocus(int32_t b, int32_t e) : beg1(b), end0(e) {}

  // compare between genomeLocus
  // l : rhs argument of the same type
  // returns true iff this < l
  inline bool operator< (const posLocus& l) const {
    if ( beg1 == l.beg1 ) { return end0 < l.end0; }
    else return beg1 < l.beg1;    
  }

  // compare between genomeLocus
  // l : rhs argument of the same type
  // returns true iff this == l
  inline bool operator== (const posLocus& l) const {
    return ( ( beg1 == l.beg1 ) && ( end0 == l.end0 ) );
  }
  
  // returns length of the interval
  int32_t length() const { return end0-beg1+1; }

  // returns the total overlapping intervals
  int32_t overlapBases(int32_t _beg1, int32_t _end0) const {
    if ( ( beg1 <= end0 ) && ( _beg1 <= end0 ) ) {
      return ( _beg1 < end0 ? _beg1 : end0 ) - ( beg1 < end0 ? end0 : beg1 ) + 1;
    }
    else return 0;    
  }

  int32_t overlapBases(const posLocus& l) const {
    return overlapBases(l.beg1, l.end0);
  }

  bool overlaps(int32_t _beg1, int32_t _end0) const {
    if ( ( beg1 <= _end0 )  && ( _beg1 <= end0 ) )
      return true;
    else 
      return false;
  }

  // check overlap with other locus
  bool overlaps (const posLocus& l) const {
    return overlaps(l.beg1, l.end0);
  }

  // merge two locus if possible
  bool merge (const posLocus& l) {
    if ( ( beg1-1 <= l.end0 )  && ( l.beg1-1 <= end0 ) ) {
      if ( l.beg1 < beg1 ) beg1 = l.beg1;
      if ( l.end0 > end0 ) end0 = l.end0;
      return true;
    }
    else {
      return false;
    }
  }

  // check whether the interval contains a particular position in 1-based coordinate
  bool contains1(int32_t pos1 = INT_MAX) const {
    return ( ( pos1 >= beg1 ) && ( pos1 <= end0 ) );
  }

  // check whether the interval contains a particular position in 0-based coordinate  
  bool contains0(int32_t pos0) const { return contains1(pos0+1); }

  // parse a string in [chr]:[beg1]-[end0] format 
  static bool parseRegion(const char* region, std::string& chrom, int32_t& beg1, int32_t& end0) {
    char buf[255];
    strcpy(buf,region);
    const char* pcolon = strchr(region,':');
    const char* pminus = strchr(pcolon+1,'-');
    if ( pcolon == NULL )
      error("Cannot parse %s in genomeLocus::genomeLocus()");
    chrom.assign(region,0,pcolon-region);
    beg1 = atoi(pcolon+1);
    if ( pminus == NULL ) end0 = INT_MAX;
    else {
      end0 = atoi(pminus+1);
      if ( end0 == 0 ) end0 = INT_MAX;
    }
    return true;
  }
};

class gtfGene;
class gtfTranscript;
class gtfCDS;
class gtfElement;
class gtfExonElement;

class gtfElement {
public:
  gtfElement* parent;
  std::string type;
  posLocus    locus;
  gtfElement(int32_t _start, int32_t _end, const char* _type, gtfElement* _parent);

  //virtual void printElement();
};

struct gtfComp {
  bool operator()(const gtfElement* lhs, const gtfElement* rhs) const {
    if ( lhs->locus == rhs->locus ) {
      return ((int64_t)rhs - (int64_t)lhs > 0);
    }
    else return ( lhs->locus < rhs->locus );
  }
};

class gtfCDS : public gtfElement {
public:
  uint8_t frame;
  gtfCDS(int32_t _start, int32_t _end, const char* sframe, gtfElement* _parent);

  //virtual void printElement();  
};

class gtfGene : public gtfElement {
public:
  std::string geneId;
  std::string geneName;
  std::string geneType;
  std::set<gtfTranscript*,gtfComp> transcripts;
  std::string seqname;
  bool fwdStrand;

  gtfGene(const char* _seqname, int32_t _start, int32_t _end, bool _fwdStrand, std::string& _gid, std::string& _gname, std::string& _gtype);

  //virtual void printElement();  
};

class gtfTranscript : public gtfElement {
public:
  std::string transcriptId;
  std::string transcriptType;
  std::set<gtfElement*,gtfComp> exons;
  std::set<gtfCDS*,    gtfComp> CDSs;
  std::set<gtfElement*,gtfComp> UTRs;
  std::set<gtfElement*,gtfComp> start_codons;
  std::set<gtfElement*,gtfComp> stop_codons;    

  gtfTranscript(int32_t _start, int32_t _end, std::string& _tid, std::string& _ttype, gtfElement* _parent);

  //virtual void printElement();  
};

// An object that represent a GTF file
class gtf {
// transcript : gene
public:
  gtf(const char* gtfFile, bool proteinCodingOnly = false, bool addChrPrefix = false, bool removeChrPrefix = false);
  ~gtf();

  int32_t maxGeneLength;
  int32_t maxTranscriptLength;
  int32_t maxExonLength;
  int32_t maxCDSLength;
  int32_t maxUTRLength;
  int32_t maxStartCodonLength;
  int32_t maxStopCodonLength;

  typedef std::multimap<posLocus, gtfElement*> gtf_chr_t;
  typedef std::map<std::string, gtf_chr_t >::iterator gtf_chr_it_t;
  std::map<std::string, gtf_chr_t > mmap;

  typedef gtfIntervalTree<int32_t,gtfElement*> gtf_ivt_t;
  //typedef gtfInterval<int32_t,gtfElement*> gtf_iv_t;  
  std::map<std::string, gtf_ivt_t> chr2ivt;
  
  std::map<std::string, gtfGene*> gid2Gene;
  std::map<std::string, gtfTranscript*> tid2Transcript;
  std::map<std::string, gtfGene*> tid2Gene;

  typedef std::multimap<posLocus, gtfElement*>::iterator gtf_elem_it_t;  

  gtf_chr_it_t   curChrIt;
  gtf_elem_it_t  curElemIt;

  void rewind();
  bool next();
  bool isend() const;

  bool addGene(const char* seqname, int32_t start, int32_t end,
	       const char* strand, 
	       std::string& gid, std::string& gname, std::string& gtype);

  bool addTranscript(const char* seqname, int32_t start, int32_t end,
	       const char* strand, 
	       std::string& gid, std::string& tid, std::string& ttype);

  bool addExon(const char* seqname, int32_t start, int32_t end,
	       const char* strand,  
	       std::string& tid);

  bool addUTR(const char* seqname, int32_t start, int32_t end,
	      const char* strand, 
	       std::string& tid);

  bool addCDS(const char* seqname, int32_t start, int32_t end,
	       const char* strand, const char* frame,
	       std::string& tid);    

  bool addStartCodon(const char* seqname, int32_t start, int32_t end,
	       const char* strand,
	       std::string& tid);

  bool addStopCodon(const char* seqname, int32_t start, int32_t end,
	       const char* strand, 
	       std::string& tid);

  bool checkTranscriptSanity(std::string& tid, const char* seqname, const char* strand);

  int32_t findOverlappingElements(const char* seqname, int32_t start, int32_t end, std::set<gtfElement*>& results);
};

#endif
