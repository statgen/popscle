#ifndef __GENOME_LOCI_H
#define __GENOME_LOCI_H

#include <vector>
#include <set>
#include <map>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <algorithm>

#include "Error.h"
#include "hts_utils.h"

// A single genomic int32_terval
class genomeLocus {
 public:
  std::string chrom; // chromosome name
  int32_t beg1; // includes 1-based, excludes 0-based
  int32_t end0; // excludes 1-based, includes 0-based
  char buf[255];

  genomeLocus(const char* c, int32_t b, int32_t e) : chrom(c), beg1(b), end0(e) {
    sprintf(buf,"%s:%d-%d",c,b,e);
  }

  // convert [chr]:[beg1]-[end0] string int32_to int32_terval
  // 20:100-110 means [100,110] in 1-based [100,111) in 1-based [99,110) in 0-based
  genomeLocus(const char* region) {
    strcpy(buf,region);
    const char* pcolon = strchr(region,':');
    const char* pminus = strchr(pcolon+1,'-');
    //if ( ( pcolon == NULL ) || ( pminus == NULL ) )
    if ( pcolon == NULL )
      error("Cannot parse %s in genomeLocus::genomeLocus()");
    chrom = std::string(region,0,pcolon-region);
    beg1 = atoi(pcolon+1);
    if ( pminus == NULL ) end0 = INT_MAX;
    else {
      end0 = atoi(pminus+1);
      if ( end0 == 0 ) end0 = INT_MAX;
    }
    //notice("%s %s %d %d", region, chrom.c_str(), beg1, end0);
  }

  const char* toString() const { 
    return buf;
  }

  // compare between genomeLocus
  bool operator< (const genomeLocus& l) const {
    if ( chrom == l.chrom ) {
      if ( beg1 == l.beg1 ) {
	return ( end0 < l.end0 );
      }
      else {
	return ( beg1 < l.beg1 );
      }
    }
    else {
      int32_t n1 = atoi(chrom.c_str());
      int32_t n2 = atoi(l.chrom.c_str());
      if ( ( n1 == 0 ) && ( n2 == 0 ) ) {
	return chrom < l.chrom;
      }
      else if ( ( n1 > 0 ) && ( n2 > 0 ) ) {
	return n1 < n2;
      }
      else { // treat n1 == 0 as infinite
	return ( n1 > 0 ) ? true : false;
      }
    }
  }

  // length
  unsigned long length() const { return end0-beg1+1; }

  int32_t overlapBases(const char* _chrom, int32_t _beg1, int32_t _end0) const {
    if ( chrom == _chrom ) {
      if ( ( beg1 <= _end0 ) && ( _beg1 <= end0 ) ) {
	return ( _beg1 < end0 ? _beg1 : end0 ) - ( beg1 < end0 ? end0 : beg1 ) + 1;
      }
      else {
	return 0;
      }
    }
    else {
      return 0;
    }    
  } 

  bool overlaps(const char* _chrom, int32_t _beg1, int32_t _end0) const {
    if ( !chrom.empty() && ( chrom == _chrom ) ) {
      if ( ( beg1 <= _end0 )  && ( _beg1 <= end0 ) ) {
	return true;
      }
      else {
	return false;
      }
    }
    else {
      return false;
    }    
  }

  // check overlap with other locus
  bool overlaps (const genomeLocus& l) const {
    if ( chrom == l.chrom ) {
      if ( ( beg1 <= l.end0 )  && ( l.beg1 <= end0 ) ) {
	return true;
      }
      else {
	return false;
      }
    }
    else {
      return false;
    }
  }

  // merge two locus if possible
  bool merge (const genomeLocus& l) {
    if ( chrom == l.chrom ) {
      if ( ( beg1-1 <= l.end0 )  && ( l.beg1-1 <= end0 ) ) {
	if ( l.beg1 < beg1 ) beg1 = l.beg1;
	if ( l.end0 > end0 ) end0 = l.end0;
	return true;
      }
      else {
	return false;
      }
    }
    else {
      return false;
    }
  }

  // check if it contains pos
  bool contains0(const char* chr, int32_t pos0) const { return contains1(chr,pos0+1); }

  bool contains1(const char* chr = NULL, int32_t pos1 = INT_MAX) const {
    notice("contains1 called pos1 = %d", pos1);
    notice("chr=%s",chr);
    if ( ( chr == NULL ) || ( chrom == chr ) ) {
      return ( ( pos1 >= beg1 ) && ( pos1 <= end0 ) );
    }
    else {
      return false;
    }
  }
};

// Collection of genomic locus
class genomeLoci {
 public:
  std::set<std::string> chroms;
  std::set<genomeLocus> loci;
  std::set<genomeLocus>::iterator it;
  bool overlapResolved;
  int32_t maxLength;

  genomeLoci() : overlapResolved(false), maxLength(0) {}
  genomeLoci(const char* reg) : overlapResolved(false), maxLength(0) {
    add(reg);
    resolveOverlaps();
  }

  // functions for iterating each locus
  inline void rewind() { it = loci.begin(); }
  inline bool next() { ++it; return ( it != loci.end() ); }
  inline bool isend() { return ( it == loci.end() );  }
  inline const genomeLocus& currentLocus() { return (*it); }

  // check the size 
  bool empty() { return loci.empty(); }
  
  bool clear() {
    if ( loci.empty() ) return false;
    
    loci.clear();
    it = loci.begin();
    overlapResolved = false;
    maxLength = 0;
    return true;
  }
 
  int32_t numLocus() const { return (int32_t)loci.size(); }

  bool openBED(const char* file) {
    clear();
    
    htsFile* fp = hts_open(file, "r");
    if ( fp == NULL )
      error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, file);
    
    kstring_t str = {0,0,0};
    int32_t lstr = 0;
    int32_t nfields = 0;
    int32_t* fields = NULL;
    int32_t i;
    // model list is assumed to have [INFO_KEY] [MODEL_FILE] [INFO_DESCRIPTION = INFO_KEY if empty]
    for( i=0; ( lstr = hts_getline(fp, KS_SEP_LINE, &str) ) >= 0; ++i ) {
      if ( fields != NULL ) { free(fields); fields = NULL; } // free the fields once allocated
      fields = ksplit(&str, 0, &nfields);
      if ( nfields < 3 )
	error("[E:%s:%d %s] Less than three columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, file);      
      // typically bed files have beg0 and end0
      add(&str.s[fields[0]], atoi(&str.s[fields[1]])+1, atoi(&str.s[fields[2]]));
      //if ( i % verbose == 0 )
      //notice("Processing %d lines to mask at %s:%d-%d in %s..", i, &str.s[fields[0]], &str.s[fields[1]], &str.s[fields[2]], file);
    }
    hts_close(fp);
    
    notice("Processed %d lines from %s, maxLength = %d", i, file, maxLength);
    
    resolveOverlaps();

    notice("After removing overlaps, %d intervals remained, maxLength = %d, total length = %u", (int32_t)loci.size(), maxLength, totalLength());

    rewind();

    return i > 0;
  }

  // add a locus
  bool add(const char* chr, int32_t beg1, int32_t end0) {
    overlapResolved = false;
    if ( end0-beg1+1 > maxLength ) maxLength = end0-beg1+1;
    std::pair<std::set<genomeLocus>::iterator, bool> ret = loci.insert(genomeLocus(chr,beg1,end0));
    it = ret.first;
    if ( ret.second )
      chroms.insert(ret.first->chrom);

    //notice("bar %s", loci.begin()->toString());
    return ret.second;
  }

  // add a locus
  bool add(const char* region) {
    overlapResolved = false;
    std::pair<std::set<genomeLocus>::iterator, bool> ret = loci.insert(genomeLocus(region));
    if ( ret.second )
      chroms.insert(ret.first->chrom);
    it = ret.first;
    int32_t l = ret.first->end0 - ret.first->beg1 + 1;
    if ( l > maxLength ) maxLength = l;

    //notice("bar %s", loci.begin()->toString());
    
    return ret.second;
  }

  // Resolve overlapping int32_tervals
  int32_t resolveOverlaps() {
    if ( !overlapResolved ) {
      std::set<genomeLocus>::iterator it;
      std::set<genomeLocus>::iterator prev;
      int32_t numMerged = 0;
      for(it = loci.begin(); it != loci.end(); ++it) {
	if ( it != loci.begin() ) {
	  if ( prev->overlaps(*it) ) {
	    // if overlaps, erase both and insert merged one
	    genomeLocus locus = *prev;
	    locus.merge(*it);
	    if ( (int32_t)locus.length() > maxLength ) maxLength = locus.length();
	    loci.erase(it);
	    loci.erase(prev);
	    prev = it = loci.insert(locus).first;
	    ++numMerged;
	  }
	  else {
	    prev = it;
	  }
	}
	else {
	  prev = it;
	}
      }
      overlapResolved = true;
      return numMerged;
    }
    else {
      return 0;
    }
    return 0;
  }

  unsigned long totalLength() const {
    //resolveOverlaps();
    unsigned long sz = 0;
    std::set<genomeLocus>::iterator it2;
    for(it2 = loci.begin(); it2 != loci.end(); ++it2) {
      sz += it2->length();
    }
    return sz;
  }

  bool hasChrom(const char* chr) {
    return ( chroms.find(chr) != chroms.end() );
  }

  bool moveTo(const char* chr = NULL, int32_t pos1 = INT_MAX) {
    notice("[%s:%d %s] (%s, %d)", __FILE__, __LINE__, __PRETTY_FUNCTION__, chr == NULL ? "NULL" : chr, pos1);
    
    if ( loci.empty() ) return false;

    if ( ( it != loci.end() ) && it->contains1(chr, pos1) ) return true;
    
    if ( chr == NULL ) chr = it->chrom.c_str();
    
    genomeLocus locus(chr, pos1, pos1);
    it = loci.lower_bound(locus);
    if ( it == loci.begin() ) { // do nothing
      notice("beg");
      return (it->contains1(chr,pos1));
    }
    else if ( it == loci.end() ) {
      notice("end");      
      std::set<genomeLocus>::iterator i = it;
      --i;
      if ( i->contains1(chr,pos1) ) { it = i; return true; }
      else { rewind(); return false; }
    }
    else {
      notice("mid");                  
      if ( it->contains1(chr,pos1) ) return true;
      else {
	std::set<genomeLocus>::iterator i = it;
	--i;
	if ( i->contains1(chr,pos1) ) { it = i; return true; }
	else { rewind(); return false; }
      }
    }
  }

  bool contains1(const char* chr, int32_t pos1) {
    if ( loci.empty() ) return false;
    notice("contains1(%s,%d) called", chr, pos1);    
    genomeLocus locus(chr, pos1, pos1);
    std::set<genomeLocus>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->chrom != chr ) ++it2;    
    while( it2 != loci.end() && ( it2->chrom == chr ) && ( it2->beg1 <= pos1 ) ) {
      if ( it2->end0 >= pos1 ) return true;
      ++it2;
    }
    return false;
  }

  bool overlaps(const char* chr, int32_t beg1, int32_t end0) {
    if ( loci.empty() ) return false;
    
    genomeLocus locus(chr, overlapResolved ? beg1 : beg1-maxLength, overlapResolved ? beg1 : beg1-maxLength);
    if ( loci.empty() ) return false;
    std::set<genomeLocus>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->chrom != chr ) ++it2;
    while( it2 != loci.end() && ( it2->chrom == chr ) && ( it2->beg1 <= end0 ) ) {
      if ( ( it2->beg1 <= end0 ) && ( beg1 <= it2->end0 ) )
	return true;
      ++it2;
    }
    //notice("%s:%d-%d",it2->chrom.c_str(),it2->beg1,it2->end0);
    return false;
  }

  bool contains(const char* chr, int32_t beg1, int32_t end0) {
    if ( loci.empty() ) return false;

    resolveOverlaps();
    genomeLocus locus(chr, beg1-maxLength, beg1-maxLength);
    std::set<genomeLocus>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->chrom != chr ) ++it2;    
    while( it2 != loci.end() && ( it2->chrom == chr ) && ( it2->beg1 <= end0 ) ) {
      if ( ( it2->beg1 <= beg1 ) && ( end0 <= it2->end0 ) )
	return true;
      ++it2;
    }
    return false;    
  }
};

// Collection of genomic locus
template <class T>
class genomeLocusMap {
 public:
  std::set<std::string> chroms;  
  std::map<genomeLocus,T> loci;
  typename std::map<genomeLocus,T>::iterator it;
  int32_t maxLength;

 genomeLocusMap() : maxLength(0) { it = loci.end(); }
  genomeLocusMap(const char* reg, const T& val) : maxLength(0) {
    add(reg, val);
    it = loci.end();
  }

  // functions for iterating each locus
  void rewind() { it = loci.begin(); }
  bool next() { ++it; return ( it != loci.end() ); }
  bool isend() { return ( it == loci.end() ); }
  const genomeLocus& currentLocus() { return (*it); }

  // check the size 
  bool empty() { return loci.empty(); }
  
  bool clear() {
    if ( loci.empty() ) return false;
    loci.clear();
    it = loci.begin();
    maxLength = 0;
    return true;
  }
    
  int32_t numLocus() const { return (int32_t)loci.size(); }

  // add a locus
  bool add(const char* chr, int32_t beg1, int32_t end0, const T& val) {
    if ( end0-beg1+1 > maxLength ) maxLength = end0-beg1+1;
    std::pair<typename std::map<genomeLocus,T>::iterator, bool> ret = loci.insert(std::pair<genomeLocus,T>(genomeLocus(chr,beg1,end0),val));
    it = ret.first;
    if ( ret.second )
      chroms.insert(chr);
    return ret.second;
  }
  
  // add a locus
  bool add(const char* region, const T& val) {
    std::pair<typename std::map<genomeLocus,T>::iterator, bool> ret = loci.insert(std::pair<genomeLocus,T>(region,val));
    int32_t l = ret.first->end0 - ret.first->beg1 + 1;
    if ( ret.second )
      chroms.insert(ret.first->chrom);    
    if ( l > maxLength ) maxLength = l;
    return ret.second;
  }

  unsigned long totalLength() const {
    unsigned long sz = 0;
    typename std::map<genomeLocus,T>::iterator it2;
    for(it2 = loci.begin(); it2 != loci.end(); ++it2) {
      sz += it2->first.length();
    }
    return sz;
  }

  bool moveTo(const char* chr = NULL, int32_t pos1 = INT_MAX) {
    notice("[%s:%d %s] (%s, %d)", __FILE__, __LINE__, __PRETTY_FUNCTION__, chr == NULL ? "NULL" : chr, pos1);

    if ( loci.empty() ) return false;    

    if ( ( it != loci.end() ) && it->first.contains1(chr, pos1) ) return true;    
    
    if ( chr == NULL ) chr = it->first.chrom.c_str();

    genomeLocus locus(chr, pos1, pos1);
    it = loci.lower_bound(locus);
    if ( it == loci.begin() ) { // do nothing
      return (it->first.contains1(chr,pos1));
    }
    else if ( it == loci.end() ) {
      typename std::map<genomeLocus,T>::iterator i = it;
      --i;
      if ( i->first.contains1(chr,pos1) ) { it = i; return true; }
      else { rewind(); return false; }
    }
    else {
      if ( it->first.contains1(chr,pos1) ) return true;
      else {
	typename std::map<genomeLocus,T>::iterator i = it;
	--i;
	if ( i->first.contains1(chr,pos1) ) { it = i; return true; }
	else { rewind(); return false; }
      }
    }
  }

  bool contains1(const char* chr, int32_t pos1) {
    notice("contains1(%s,%d) called", chr, pos1);
    genomeLocus locus(chr, pos1, pos1);
    typename std::map<genomeLocus,T>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->chrom != chr ) ++it2;    
    while( it2 != loci.end() && ( it2->first.chrom == chr ) && ( it2->first.beg1 <= pos1 ) ) {
      if ( it2->first.end0 >= pos1 ) return true;
      ++it2;
    }
    return false;
  }

  bool overlaps(const char* chr, int32_t beg1, int32_t end0) {
    genomeLocus locus(chr, beg1-maxLength, beg1-maxLength);
    if ( loci.empty() ) return false;
    typename std::map<genomeLocus,T>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->first.chrom != chr ) ++it2;
    while( it2 != loci.end() && ( it2->first.chrom == chr ) && ( it2->first.beg1 <= end0 ) ) {
      if ( ( it2->first.beg1 <= end0 ) && ( beg1 <= it2->first.end0 ) )
	return true;
      ++it2;
    }
    return false;
  }

  bool contains(const char* chr, int32_t beg1, int32_t end0) {
    if ( loci.empty() ) return false;

    genomeLocus locus(chr, beg1-maxLength, beg1-maxLength);
    typename std::map<genomeLocus,T>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->first.chrom != chr ) ++it2;    
    while( it2 != loci.end() && ( it2->first.chrom == chr ) && ( it2->first.beg1 <= end0 ) ) {
      if ( ( it2->first.beg1 <= beg1 ) && ( end0 <= it2->first.end0 ) )
	return true;
      ++it2;
    }
    return false;    
  }
};

#endif
