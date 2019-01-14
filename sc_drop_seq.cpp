#include <cmath>
#include "sc_drop_seq.h"
#include "PhredHelper.h"

double logAdd(double la, double lb) {
  if ( la > lb ) { return la + log(1.0 + exp(lb-la)); }
  else           { return lb + log(1.0 + exp(la-lb)); }  
}

int32_t sc_dropseq_lib_t::add_snp(int32_t _rid, int32_t _pos, char _ref, char _alt, double _af, double* _gps) {
  snps.resize(nsnps+1);
  snp_umis.resize(nsnps+1);
  sc_snp_t& snp = snps.back();
  snp.rid = _rid;
  snp.pos = _pos;
  snp.ref = _ref;
  snp.alt = _alt;
  snp.af  = _af;
  snp.gps = _gps;

  //if ( _pos % 1000 == 0 ) notice("%lf %lf %lf %lf %lf %lf",_gps[3],_gps[4],_gps[5],_gps[6],_gps[7],_gps[8]);
  
  ++nsnps;
  return nsnps-1;
}

int32_t sc_dropseq_lib_t::add_cell(const char* barcode) {
  std::map<std::string,int32_t>::iterator it = bc_map.find(barcode);
  if ( it == bc_map.end() ) {
    bc_map[barcode] = nbcs;
    bcs.push_back(barcode);
    
    cell_umis.resize( nbcs + 1 );
    cell_totl_reads.resize( nbcs + 1 );
    cell_pass_reads.resize( nbcs + 1 );
    cell_uniq_reads.resize( nbcs + 1 );
    cell_scores.resize(nbcs + 1);
    ++nbcs;
    
    return (nbcs-1);
  }
  else return it->second;
}

bool sc_dropseq_lib_t::add_read(int32_t snpid, int32_t cellid, const char* umi, char allele, char qual) {
  std::map<int32_t,sc_snp_droplet_t*>::iterator it = snp_umis[snpid].find(cellid);
  sc_snp_droplet_t* p_snp_drop = NULL;
  bool ret = false;

  //if ( ( allele > 2 ) || ( qual > 50 ) )
  //  error("umi = %s, allele = %c, qual = %c", umi, allele, qual);

  ++cell_pass_reads[cellid];

  // check if (snp,cell) is empty
  if ( it == snp_umis[snpid].end() ) {
    p_snp_drop = new sc_snp_droplet_t;
    (*p_snp_drop)[umi] = (uint32_t)( ( allele << 24 ) | ( qual << 16 ) | 0x01 );
    snp_umis[snpid][cellid] = p_snp_drop;
    //p_snp_drop = snp_umis[snpid][cellid];
    ret = true;
  }
  // check if (snp,cell,umi) is empty
  else {
    sc_snp_droplet_it_t it2 = it->second->find(umi);
    if ( it2 == it->second->end() ) {
      (*it->second)[umi] = (uint32_t)( ( allele << 24 ) | ( qual << 16 ) | 0x01 );
      ret = true;
    }
    else {
      ++(it->second->at(umi));
    }
    p_snp_drop = it->second;
  }

  std::map<int32_t,sc_snp_droplet_t*>::iterator it3 = cell_umis[cellid].find(snpid);
  if ( it3 == cell_umis[cellid].end() ) {
    cell_umis[cellid][snpid] = p_snp_drop;
  }
  else if ( it3->second != p_snp_drop ) {
    notice("[E:%s] Conflict : Multiple (cellid,snpid) pair %d %d, %x %x", __PRETTY_FUNCTION__, cellid, snpid, it3->second, p_snp_drop);
    notice("%s",p_snp_drop->begin()->first.c_str());
    notice("%d",it3->first);    
    notice("%x",it3->second);
    notice("%u",it3->second->size());
    notice("%s",it3->second->begin()->first.c_str());        
    error("[E:%s] Conflict : Multiple (cellid,snpid) pair %d %d, %s %s", __PRETTY_FUNCTION__, cellid, snpid, p_snp_drop->begin()->first.c_str(), it3->second->begin()->first.c_str());
  }
  if ( ret ) ++cell_uniq_reads[cellid];
  return ret;
}


int32_t sc_dropseq_lib_t::load_from_plp(const char* plpPrefix, BCFFilteredReader* pvr, const char* field, double genoErrorOffset, double genoErrorCoeffR2, const char* r2info, bool loadUMI) {
  if ( loadUMI == true )
    error("[E:%s] loadUMI = true option is not implemented yet", __PRETTY_FUNCTION__);

  int32_t nv = 0;
  if ( pvr != NULL ) { // variant sites are provided..
    if ( pvr->read() == NULL )  // attempt to read a variant from VCF
      error("[E:%s Cannot read any single variant from %s]", __PRETTY_FUNCTION__, pvr->bcf_file_name.c_str());

    // attempt to parse genotype field to make sure the field exists
    if ( ! pvr->parse_posteriors(pvr->cdr.hdr, pvr->cursor(), field, 0) )
      error("[E:%s] Cannot parse posterior probability at %s:%d", __PRETTY_FUNCTION__, bcf_hdr_id2name(pvr->cdr.hdr,pvr->cursor()->rid), pvr->cursor()->pos+1);

    nv = pvr->get_nsamples();    
  }

  verbose(70, "Loading pileup information with prefix %s", plpPrefix);
  verbose(50, "Reading barcode information from %s.cel.gz..", plpPrefix);  
  char fname[65535];

  // reading the droplet information first..
  sprintf(fname, "%s.cel.gz", plpPrefix);
  tsv_reader tsv_bcdf(fname);

  std::vector<int32_t> tmp_cell_totl_reads;
  std::vector<int32_t> tmp_cell_totl_umis;  
  std::vector<int32_t> tmp_cell_uniq_reads;
  std::vector<int32_t> tmp_cell_num_snps;    

  int32_t n_expected_toks = 6;
  if ( tsv_bcdf.read_line() > 0 ) {
    if ( ( tsv_bcdf.nfields == 5 ) &&   // for backward compatibility
	 ( ( strcmp("#DROPLET_ID",tsv_bcdf.str_field_at(0)) != 0 ) ||
	   ( strcmp("BARCODE",tsv_bcdf.str_field_at(1)) != 0 ) ||
	   ( strcmp("NUM.READ",tsv_bcdf.str_field_at(2)) != 0 ) ||
	   ( strcmp("NUM.UMI",tsv_bcdf.str_field_at(3)) != 0 ) ||
	   ( strcmp("NUM.SNP",tsv_bcdf.str_field_at(4)) != 0 ) ) ) {
      error("The header line of %s.cel.gz is malformed or outdated. Expecting #DROPLET_ID BARCODE NUM.READ NUM.UMI NUM.SNP", plpPrefix);
    }
    else if ( ( tsv_bcdf.nfields == 6 ) &&
	 ( ( strcmp("#DROPLET_ID",tsv_bcdf.str_field_at(0)) != 0 ) ||
	   ( strcmp("BARCODE",tsv_bcdf.str_field_at(1)) != 0 ) ||
	   ( strcmp("NUM.READ",tsv_bcdf.str_field_at(2)) != 0 ) ||
	   ( strcmp("NUM.UMI",tsv_bcdf.str_field_at(3)) != 0 ) ||
	   ( strcmp("NUM.UMIwSNP",tsv_bcdf.str_field_at(4)) != 0 ) ||	   
	   ( strcmp("NUM.SNP",tsv_bcdf.str_field_at(5)) != 0 ) ) ) {
      error("The header line of %s.cel.gz is malformed or outdated. Expecting #DROPLET_ID BARCODE NUM.READ NUM.UMI NUM.UMIwSNP NUM.SNP", plpPrefix);
    }
    else if ( ( tsv_bcdf.nfields < 5 ) || ( tsv_bcdf.nfields > 6 ) ) {
      error("The header line of %s.cel.gz is malformed or outdated. Expecting #DROPLET_ID BARCODE NUM.READ NUM.UMI (NUM.UMIwSNP-optional) NUM.SNP", plpPrefix);  
    }
    n_expected_toks = tsv_bcdf.nfields;
  }
  else error("Cannot read the first line of %s.cel.gz", plpPrefix);

  while( tsv_bcdf.read_line() > 0 ) {
    int32_t new_id = add_cell(tsv_bcdf.str_field_at(1));
    if ( new_id != tsv_bcdf.int_field_at(0) )
      error("[E:%s] Observed DROPLET_ID %d is different from expected DROPLET_ID. Did you modify the digital pileup files by yourself?", __PRETTY_FUNCTION__, tsv_bcdf.int_field_at(0), new_id);
    tmp_cell_totl_reads.push_back(tsv_bcdf.int_field_at(2));
    if ( n_expected_toks == 5 ) {
      tmp_cell_uniq_reads.push_back(tsv_bcdf.int_field_at(3));
      tmp_cell_num_snps.push_back(tsv_bcdf.int_field_at(4));
    }
    else {
      tmp_cell_totl_umis.push_back(tsv_bcdf.int_field_at(3));  // could be zero. no sanity check on this.      
      tmp_cell_uniq_reads.push_back(tsv_bcdf.int_field_at(4));
      tmp_cell_num_snps.push_back(tsv_bcdf.int_field_at(5));      
    }
  }
  verbose(50, "Finished loading %d droplets..", nbcs);

  // reading the variant information next..
  verbose(50, "Reading variant information from %s.var.gz..", plpPrefix);  
  sprintf(fname, "%s.var.gz", plpPrefix);
  tsv_reader tsv_varf(fname);  
  if ( tsv_varf.read_line() > 0 ) {
    if ( ( tsv_varf.nfields != 6 ) ||
	 ( strcmp("#SNP_ID",tsv_varf.str_field_at(0)) != 0 ) ||
	 ( strcmp("CHROM",tsv_varf.str_field_at(1)) != 0 ) ||
	 ( strcmp("POS",tsv_varf.str_field_at(2)) != 0 ) ||
	 ( strcmp("REF",tsv_varf.str_field_at(3)) != 0 ) ||
	 ( strcmp("ALT",tsv_varf.str_field_at(4)) != 0 ) ||	   
	 ( strcmp("AF",tsv_varf.str_field_at(5)) != 0 ) )
      error("THe header line of %s.var.gz is malformed or outdated. Expecting #SNP_ID CHROM POS REF ALT AF", plpPrefix);
  }
  else error("Cannot read the first line of %s.var.gz", plpPrefix);

  int nrd = 0;

  float* r2flts = NULL; // for extracting R2 fields..
  int32_t n_r2flts = 0;
  
  while( tsv_varf.read_line() > 0 ) {
    const char* chr = tsv_varf.str_field_at(1);
    if ( chr2rid.find(chr) == chr2rid.end() ) {
      int32_t newrid = chr2rid.size();
      chr2rid[chr] = newrid;
      rid2chr.push_back(chr);
    }
    int32_t rid = chr2rid[chr];
    int32_t pos = tsv_varf.int_field_at(2);
    char    ref = tsv_varf.str_field_at(3)[0];
    char    alt = tsv_varf.str_field_at(4)[0];
    double  af  = tsv_varf.double_field_at(5);

    if ( pvr == NULL ) { // no VCFs were provided as argument
      if ( add_snp(rid, pos, ref, alt, af, NULL) + 2 != tsv_varf.nlines )
	error("Expected SNP nID = %d but observed %d", tsv_varf.nlines-1, nsnps-1);
    }
    else {
      // find the variant from VCF. The VCF must be overlapping with what was provided before
      bcf1_t* v = pvr->cursor(); // read it and have not iterated over yet

      if ( rand() % 10000 == 0 )
	verbose(50, "Reading variant info %s:%d:%c:%c at %s:%d:%c:%c", tsv_varf.str_field_at(0), pos, ref, alt, bcf_hdr_id2name(pvr->cdr.hdr, v->rid), v->pos+1, v->d.allele[0][0], v->d.allele[1][0] );
      
      // we need to
      // 1. If the position exactly matches
      //    1.1.  If the alleles do not match, throw error
      //    1.2.  Otherwise, add the cursor (GP/GT should be available)
      // 2. If the cursor already passed the current position from var.gz
      //    Pretend that there is a SNP, but with empty GP (add a way to handle empty GP)
      // 3. If the cursor has not passed the current position from var.gz
      //    Simply ignore the and iterate more markers from VCF
      bool found = false;
      bool passed = false;
      while ( ! ( found || passed ) ) {
	if ( pvr->eof ) { passed = true; }
	else if ( v->rid > rid ) { passed = true; }
	else if ( v->rid == rid ) { // same chromosome
	  if ( v->pos + 1 > pos ) { passed = true; }
	  else if ( v->pos + 1 == pos ) {
	    if ( ( v->d.allele[0][0] != ref ) || ( v->d.allele[1][0] != alt ) )
	      //error("Could not find variant %s:%d:%c:%c from genotype VCF file, cursor:%d:%d:%c:%c", tsv_varf.str_field_at(0), pos, ref, alt, v->rid, v->pos +1, v->d.allele[0][0], v->d.allele[1][0]);	  
	      //multi-allele,cursor passed, read next var.gz marker
	    {
		    passed = true;
		    //break;//because it has been added previously
	    }
	    else
		found = true;
	  }
	}

	if ( passed ) { // add just empty SNPs without any GPs
	  if ( add_snp(rid, pos, ref, alt, af, NULL) + 2 != tsv_varf.nlines )
	    error("Expected SNP nID = %d but observed %d", tsv_varf.nlines-1, nsnps-1);
	  break;
	}
	else if ( found ) {
	  // get GP, PL, or GT fields to convert into gps
	  if ( ! pvr->parse_posteriors(pvr->cdr.hdr, v, field, 0) )
	    error("[E:%s] Cannot parse posterior probability at %s:%d", __PRETTY_FUNCTION__, bcf_hdr_id2name(pvr->cdr.hdr,v->rid), v->pos+1);
	  double* gps = new double[nv*3];
	  double avgGPs[3] = {1e-10,1e-10,1e-10}; // represents empirical GP averages across sampleso
	  // get genotype probabilities
	  for(int32_t i=0; i < nv * 3; ++i) {
	    avgGPs[i%3] += (gps[i] = pvr->get_posterior_at(i));
	  }
	  // get average GPs to account for genoErrors
	  double sumGP = avgGPs[0] + avgGPs[1] + avgGPs[2];
	  avgGPs[0] /= sumGP;
	  avgGPs[1] /= sumGP;
	  avgGPs[2] /= sumGP;

	  // account for genotype errors as [Offset] + [1-Offset]*[1-R2]*[Coeff]
	  double err = genoErrorOffset;
	  if ( genoErrorCoeffR2 > 0 ) { // look for R2 INFO field
	    if ( ( bcf_get_info_float(pvr->cdr.hdr, v, r2info, &r2flts, &n_r2flts) < 0 ) || ( n_r2flts != 1 ) ) {
	      error("Cannot extract %s (1 float value) from INFO field at %s:%d. Cannot use --geno-error-coeff", r2info, chr, pos);
	    }
	    err += (1-genoErrorOffset) * (1-r2flts[0]) * genoErrorCoeffR2;
	  }

	  if ( err > 0.999 ) err = 0.999;
	  if ( err < 0 ) err = 0;	  
	  
	  if ( err > 0 ) { // if error is greater than zero, adjust it
	    for(int32_t i=0; i < nv * 3; ++i) {
	      gps[i] = (1-err) * gps[i] + err * avgGPs[ i % 3 ];
	    }	    
	  }

	  // Sanity check on the number of SNPs
	  if ( add_snp(rid, pos, ref, alt, af, gps) + 2 != tsv_varf.nlines )
	    error("Expected SNP nID = %d but observed %d", tsv_varf.nlines-1, nsnps-1);
	  break;
	}
	pvr->read();	
	v = pvr->cursor();
	++nrd;
      }
    }
  }
  verbose(50, "Finished loading %d variants..", nsnps);

  free(r2flts);  

  // reading pileup information finally...
  verbose(50, "Reading pileup information from %s.plp.gz..", plpPrefix);  
  sprintf(fname, "%s.plp.gz", plpPrefix);
  tsv_reader tsv_plpf(fname);
  if ( tsv_plpf.read_line() > 0 ) {
    if ( ( tsv_plpf.nfields != 4 ) ||
	 ( strcmp("#DROPLET_ID",tsv_plpf.str_field_at(0)) != 0 ) ||
	 ( strcmp("SNP_ID",tsv_plpf.str_field_at(1)) != 0 ) ||
	 ( strcmp("ALLELES",tsv_plpf.str_field_at(2)) != 0 ) ||
	 ( strcmp("BASEQS",tsv_plpf.str_field_at(3)) != 0 ) )
      error("THe header line of %s.cel.gz is malformed or outdated. Expecting #DROPLET_ID SNP_ID ALLELES BASEQS", plpPrefix);
  }
  else error("Cannot read the first line of %s.plp.gz", plpPrefix);      
  
  int32_t numi = 0;
  char buf[255];
  while( tsv_plpf.read_line() > 0 ) {
    const char* pa = tsv_plpf.str_field_at(2);
    const char* pq = tsv_plpf.str_field_at(3);
    int32_t l = (int32_t)strlen(pq);
    
    if ( (int32_t)strlen(pq) != l )
      error("Length are different between %s and %s", pa, pq);
    
    for(int32_t i=0; i < l; ++i) {
      sprintf(buf, "%x", numi++);
      ++cell_totl_reads[tsv_plpf.int_field_at(0)];    	
      add_read( tsv_plpf.int_field_at(1), tsv_plpf.int_field_at(0), buf, (char)(pa[i]-(char)'0'), (char)(pq[i]-(char)33) ); 
    }
  }
  verbose(50, "Finished loading %d UMIs in total..", numi);

  // sanity check on the observed counts
  for(int32_t i=0; i < nbcs; ++i) {
    if ( ( cell_uniq_reads[i] == tmp_cell_uniq_reads[i] ) &&
	 ( tmp_cell_num_snps[i] == (int32_t)cell_umis[i].size() ) ) {
      cell_totl_reads[i] = tmp_cell_totl_reads[i]; // overwrite
    }
  }

  return numi;
}

double calculate_snp_droplet_doublet_GL(sc_snp_droplet_t* ssd, double* gls, double alpha) {
  //double logdenom = 0;
  double tmp;
  gls[0] = gls[1] = gls[2] = gls[3] = gls[4] = gls[5] = gls[6] = gls[7] = gls[8] = 1.0;

  // iterate over each possible bases
  for(sc_snp_droplet_it_t it = ssd->begin(); it != ssd->end(); ++it) {
    uint8_t al = ( it->second >> 24 ) & 0x00ff;
    uint8_t bq = ( it->second >> 16 ) & 0x00ff;

    if ( al == 2 ) continue;

    // 0 : REF / REF -- 1        0
    // 1 : REF / HET -- 1-a/2    a/2
    // 2 : REF / ALT -- 1-a      a
    // 3 : HET / REF -- (1+a)/2  (1-a)/2
    // 4 : HET / HET -- 1/2      1/2
    // 5 : HET / ALT -- (1-a)/2  (1+a)/2
    // 6 : ALT / REF -- a        1-a
    // 7 : ALT / HET -- a/2      1-a/2
    // 8 : ALT / ALT -- 0        1

    gls[0] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 1.0           : 0.0          ) + phredConv.phred2Err[bq] / 4. );
    gls[1] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 1. - alpha/2. : alpha/2.     ) + phredConv.phred2Err[bq] / 4. );
    gls[2] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 1.0 - alpha   : alpha        ) + phredConv.phred2Err[bq] / 4. );
    gls[3] *= ( phredConv.phred2Mat[bq] * (al == 0 ? (1.+alpha)/2. : (1.-alpha)/2.) + phredConv.phred2Err[bq] / 4. );
    gls[4] *= ( phredConv.phred2Mat[bq] * (al == 0 ? .5            : .5           ) + phredConv.phred2Err[bq] / 4. );
    gls[5] *= ( phredConv.phred2Mat[bq] * (al == 0 ? (1.-alpha)/2. : (1.+alpha)/2.) + phredConv.phred2Err[bq] / 4. );
    gls[6] *= ( phredConv.phred2Mat[bq] * (al == 0 ? alpha         : 1.-alpha     ) + phredConv.phred2Err[bq] / 4. );
    gls[7] *= ( phredConv.phred2Mat[bq] * (al == 0 ? alpha/2.      : 1.-alpha/2.  ) + phredConv.phred2Err[bq] / 4. );
    gls[8] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 0.0           : 1.0          ) + phredConv.phred2Err[bq] / 4. );    

    tmp = gls[0] + gls[1] + gls[2] + gls[3] + gls[4] + gls[5] + gls[6] + gls[7] + gls[8];
    gls[0] /= tmp;
    gls[1] /= tmp;
    gls[2] /= tmp;
    gls[3] /= tmp;
    gls[4] /= tmp;
    gls[5] /= tmp;
    gls[6] /= tmp;
    gls[7] /= tmp;
    gls[8] /= tmp;        
    //logdenom += log(tmp);
  }

  for(int32_t i=0; i < 9; ++i) {
    if ( gls[i] < MIN_NORM_GL )
      gls[i] = MIN_NORM_GL;
  }
  tmp = gls[0] + gls[1] + gls[2] + gls[3] + gls[4] + gls[5] + gls[6] + gls[7] + gls[8];  

  gls[0] /= tmp;
  gls[1] /= tmp;
  gls[2] /= tmp;
  gls[3] /= tmp;
  gls[4] /= tmp;
  gls[5] /= tmp;
  gls[6] /= tmp;
  gls[7] /= tmp;
  gls[8] /= tmp;          
  //logdenom += log(tmp);

  //return logdenom;
  return 0;  
}

double calculate_snp_droplet_pileup(sc_snp_droplet_t* ssd, snp_droplet_pileup* sdp, double alpha) {
  double tmp;
  if ( sdp == NULL ) error("ERROR: NULL snp_droplet_pileup* as inp[ut");
  std::fill(sdp->gls, sdp->gls+9, 1.0);
  sdp->logdenom = 0;

  double* gls = sdp->gls;

  // iterate over each possible bases
  for(sc_snp_droplet_it_t it = ssd->begin(); it != ssd->end(); ++it) {
    uint8_t al = ( it->second >> 24 ) & 0x00ff;
    uint8_t bq = ( it->second >> 16 ) & 0x00ff;

    ++(sdp->nreads);
    
    if ( al > 1 ) continue;

    if ( al == 0 ) ++(sdp->nref);
    else if ( al == 1 ) ++(sdp->nalt);

    // 0 : REF / REF -- 1        0
    // 1 : REF / HET -- 1-a/2    a/2
    // 2 : REF / ALT -- 1-a      a
    // 3 : HET / REF -- (1+a)/2  (1-a)/2
    // 4 : HET / HET -- 1/2      1/2
    // 5 : HET / ALT -- (1-a)/2  (1+a)/2
    // 6 : ALT / REF -- a        1-a
    // 7 : ALT / HET -- a/2      1-a/2
    // 8 : ALT / ALT -- 0        1

    gls[0] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 1.0           : 0.0          ) + phredConv.phred2Err[bq] / 4. );
    gls[1] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 1. - alpha/2. : alpha/2.     ) + phredConv.phred2Err[bq] / 4. );
    gls[2] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 1.0 - alpha   : alpha        ) + phredConv.phred2Err[bq] / 4. );
    gls[3] *= ( phredConv.phred2Mat[bq] * (al == 0 ? (1.+alpha)/2. : (1.-alpha)/2.) + phredConv.phred2Err[bq] / 4. );
    gls[4] *= ( phredConv.phred2Mat[bq] * (al == 0 ? .5            : .5           ) + phredConv.phred2Err[bq] / 4. );
    gls[5] *= ( phredConv.phred2Mat[bq] * (al == 0 ? (1.-alpha)/2. : (1.+alpha)/2.) + phredConv.phred2Err[bq] / 4. );
    gls[6] *= ( phredConv.phred2Mat[bq] * (al == 0 ? alpha         : 1.-alpha     ) + phredConv.phred2Err[bq] / 4. );
    gls[7] *= ( phredConv.phred2Mat[bq] * (al == 0 ? alpha/2.      : 1.-alpha/2.  ) + phredConv.phred2Err[bq] / 4. );
    gls[8] *= ( phredConv.phred2Mat[bq] * (al == 0 ? 0.0           : 1.0          ) + phredConv.phred2Err[bq] / 4. );    

    tmp = 0;
    for(int32_t i=0; i < 9; ++i) tmp += gls[i];
    for(int32_t i=0; i < 9; ++i) gls[i] /= tmp;
    sdp->logdenom += log(tmp);
  }

  for(int32_t i=0; i < 9; ++i) {
    if ( gls[i] < MIN_NORM_GL )
      gls[i] = MIN_NORM_GL;
  }  
  tmp = 0;
  for(int32_t i=0; i < 9; ++i) tmp += sdp->gls[i];
  for(int32_t i=0; i < 9; ++i) gls[i] /= tmp;  
  
  sdp->logdenom += log(tmp);

  return sdp->logdenom;
}

double calculate_snp_droplet_GL(sc_snp_droplet_t* ssd, double* gls) {
  double logdenom = 0;
  double tmp;
  gls[0] = gls[1] = gls[2] = 1.0;
  for(sc_snp_droplet_it_t it = ssd->begin(); it != ssd->end(); ++it) {
    uint8_t al = ( it->second >> 24 ) & 0x00ff;
    uint8_t bq = ( it->second >> 16 ) & 0x00ff;

    if ( al == 2 ) continue;

    gls[0] *= ((al==0) ? phredConv.phred2Mat[bq] : phredConv.phred2Err[bq]/3.0);
    gls[1] *= (0.5 - phredConv.phred2Err[bq]/3.0);
    gls[2] *= ((al==1) ? phredConv.phred2Mat[bq] : phredConv.phred2Err[bq]/3.0);
    tmp = gls[0] + gls[1] + gls[2];
    gls[0] /= tmp;
    gls[1] /= tmp;
    gls[2] /= tmp;
    logdenom += log(tmp);
  }

  if ( gls[0] < MIN_NORM_GL ) gls[0] = MIN_NORM_GL;
  if ( gls[1] < MIN_NORM_GL ) gls[1] = MIN_NORM_GL;
  if ( gls[2] < MIN_NORM_GL ) gls[2] = MIN_NORM_GL;  
  tmp = gls[0] + gls[1] + gls[2];
  gls[0] /= tmp;
  gls[1] /= tmp;
  gls[2] /= tmp;
  logdenom += log(tmp);

  return logdenom;
  //return 0;
}
