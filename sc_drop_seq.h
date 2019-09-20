#ifndef __SC_DROP_SEQ_H
#define __SC_DROP_SEQ_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <stdint.h>

#include "Error.h"
#include "bcf_filtered_reader.h"


#define MIN_NORM_GL 1e-6

// read : read ID
// unique read : (barcode, snp, umi)
// snp  : snp ID

// read -> cell -> UMI


// 8bit  - BQ
// 8bit  - allele
// 16bit - count
typedef std::map<std::string,uint32_t> sc_snp_droplet_t;
typedef std::map<std::string,uint32_t>::iterator sc_snp_droplet_it_t;

class sc_snp_t {
 public:
  int32_t rid;
  int32_t pos;
  char ref;
  char alt;
  double af;
  double* gps;
};

double logAdd(double la, double lb);
//{
//  if ( la > lb ) { return la + log(1.0 + exp(lb-la)); }
//  else           { return lb + log(1.0 + exp(la-lb)); }  
//}

struct dropD {
  int32_t nsnps;
  int32_t nread1;
  int32_t nread2;
  double llk0;
  double llk2;

  dropD() : nsnps(0), nread1(0), nread2(0), llk0(0), llk2(0) {}
  /*
  dropD(int32_t _nsnps, double _llk0, double _llk2) :
    nsnps(_nsnps), nreads(_nreads), llk2(_llk2) {}

  void set(int32_t _nsnps, double _llk0, double _llk2) {
    nsnps = _nsnps;
    llk0 = _llk0;
    llk2 = _llk2;
  }
  */
};

struct snp_droplet_pileup {
  int32_t nreads;
  int32_t nref;
  int32_t nalt;
  double  gls[9];
  double  logdenom;

  snp_droplet_pileup() : nreads(0), nref(0), nalt(0) {
    std::fill(gls, gls+9, 1.0);    
    logdenom = 0;
  }

  void merge(const snp_droplet_pileup& other) {
    nreads += other.nreads;
    nref += other.nref;
    nalt += other.nalt;
    logdenom += other.logdenom;
    
    for(int i=0; i < 9; ++i) gls[i] *= other.gls[i];
    
    double tmp = 0;
    for(int i=0; i < 9; ++i) tmp += gls[i];
    
    logdenom += log(tmp);
    
    for(int i=0; i < 9; ++i) gls[i] /= tmp;
    
    for(int i=0; i < 9; ++i) {
      if ( gls[i] < MIN_NORM_GL ) {
	gls[i] = MIN_NORM_GL;
      }
    }
    tmp = 0;
    for(int i=0; i < 9; ++i) tmp += gls[i];
    logdenom += log(tmp);
    for(int i=0; i < 9; ++i) gls[i] /= tmp;    
  }

  void remove(const snp_droplet_pileup& other) {
    nreads -= other.nreads;
    nref -= other.nref;
    nalt -= other.nalt;
    logdenom -= other.logdenom;
    
    for(int i=0; i < 9; ++i) gls[i] /= other.gls[i];
    
    double tmp = 0;
    for(int i=0; i < 9; ++i) tmp += gls[i];
    
    logdenom += log(tmp);
    
    for(int i=0; i < 9; ++i) gls[i] /= tmp;
    
    for(int i=0; i < 9; ++i) {
      if ( gls[i] < MIN_NORM_GL ) {
	gls[i] = MIN_NORM_GL;
      }
    }
    tmp = 0;
    for(int i=0; i < 9; ++i) tmp += gls[i];
    logdenom += log(tmp);
    for(int i=0; i < 9; ++i) gls[i] /= tmp;    
  }  
};

class sc_dropseq_lib_t {
 public:
  // information to filter the input barcodes
  int32_t nbcs;
  int32_t nsnps;  
  int32_t minUMI;
  int32_t minSNP;
  int32_t minRead;
  int32_t minBQ;
  int32_t capBQ;
  
  std::set<std::string> valid_bcs;
  std::vector<int32_t>  index_bcs;
  
  // vector containing SNP & genotype info, index is snp_id
  std::map<std::string, int32_t> chr2rid;
  std::vector<std::string>       rid2chr;
  
  std::vector<sc_snp_t> snps;

  // mapper between barcode -> bcd_id  
  std::map<std::string,int32_t> bc_map;
  std::vector<std::string> bcs;
  
  // cell_umis[i]->[j] contains the map of UMIs overlapping with snp j in cell i
  std::vector< std::map<int32_t,sc_snp_droplet_t*> > cell_umis;

  // Number of pass-filtered reads and unique reads
  std::vector<int32_t> cell_totl_reads;  
  std::vector<int32_t> cell_pass_reads;
  std::vector<int32_t> cell_uniq_reads;
  std::vector<double>  cell_scores;
  
  std::vector< std::map<int32_t,sc_snp_droplet_t*> > snp_umis;
  int32_t add_snp(int32_t _rid, int32_t _pos, char _ref, char _alt, double _af, double* _gps);
  int32_t add_cell(const char* barcode);
  bool add_read(int32_t snpid, int32_t cellid, const char* umi, char allele, char qual);

  int32_t load_valid_barcodes(const char* bcdFile);
  int32_t load_from_plp(const char* plpPrefix, BCFFilteredReader* pvr = NULL, const char* field = "GP", double genoErrorOffset = 0.1, double genoErrorCoeffR2 = 0.0, const char* r2info = "R2", bool loadUMI = false);

  sc_dropseq_lib_t() : nbcs(0), nsnps(0), minUMI(0), minSNP(0), minRead(0), minBQ(1), capBQ(60) {}

  dropD calculate_droplet_clust_distance(std::map<int32_t,snp_droplet_pileup*> dropletPileup, std::map<int32_t,snp_droplet_pileup>& clustPileup);
};

double calculate_snp_droplet_GL(sc_snp_droplet_t* ssd, double* gls);
double calculate_snp_droplet_doublet_GL(sc_snp_droplet_t* ssd, double* gls, double alpha);
double calculate_snp_droplet_pileup(sc_snp_droplet_t* ssd, snp_droplet_pileup* sdp, double alpha);

struct sc_drop_comp_t {
  sc_dropseq_lib_t* pscl;
  sc_drop_comp_t(sc_dropseq_lib_t* p) : pscl(p) {}
  bool operator()(const int32_t& lhs, const int32_t& rhs) const {
    double cmp = pscl->cell_scores[lhs] - pscl->cell_scores[rhs];
    if ( cmp != 0 ) return cmp > 0;
    else return lhs > rhs;
  }
};

#endif
