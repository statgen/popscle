#include "cramore.h"
#include "bcf_filtered_reader.h"
#include "sam_filtered_reader.h"
#include "sc_drop_seq.h"

int32_t cmdCramMuxPileup(int32_t argc, char** argv) {
  SAMFilteredReader sr;
  BCFFilteredReader vr;
  std::string outPrefix;
  std::string tagGroup("CB");
  std::string tagUMI("UB");
  int32_t capBQ = 40;
  int32_t minBQ = 13;
  int32_t minTD = 0;
  sr.filt.exclude_flag = 0x0f04;
  sr.filt.minMQ = 20;  
  std::vector<std::string> smIDs;
  vr.verbose = 10000;
  sr.verbose = 1000000;  
  vr.vfilt.minMAC = 1;
  vr.vfilt.minCallRate = 0.5;
  vr.vfilt.maxAlleles = 2;  
  std::string groupList;
  std::string field("GT");
  int32_t minTotalReads = 0;
  int32_t minUniqReads = 0;
  int32_t minCoveredSNPs = 0;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&sr.sam_file_name, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("tag-group",&tagGroup, "Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB")
    LONG_STRING_PARAM("tag-UMI",&tagUMI, "Tag representing UMIs. For 10x genomiucs, use UB")

    LONG_PARAM_GROUP("Options for input VCF/BCF", NULL)
    LONG_STRING_PARAM("vcf",&vr.bcf_file_name, "Input VCF/BCF file, containing the AC and AN field")
    LONG_STRING_PARAM("field",&field,"FORMAT field to extract the genotype, likelihood, or posterior from")    
    LONG_INT_PARAM("min-mac",&vr.vfilt.minMAC, "Minimum minor allele frequency")
    LONG_DOUBLE_PARAM("min-callrate",&vr.vfilt.minCallRate, "Minimum call rate")    
    LONG_MULTI_STRING_PARAM("sm",&smIDs, "List of sample IDs to compare to (default: use all)")
    LONG_STRING_PARAM("sm-list",&vr.sample_id_list, "File containing the list of sample IDs to compare")        

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Output file prefix")
    LONG_INT_PARAM("sam-verbose",&sr.verbose, "Verbose message frequency for SAM/BAM/CRAM")
    LONG_INT_PARAM("vcf-verbose",&vr.verbose, "Verbose message frequency for VCF/BCF")
    
    LONG_PARAM_GROUP("Read filtering Options", NULL)
    LONG_INT_PARAM("cap-BQ", &capBQ, "Maximum base quality (higher BQ will be capped)")
    LONG_INT_PARAM("min-BQ", &minBQ, "Minimum base quality to consider (lower BQ will be skipped)")
    LONG_INT_PARAM("min-MQ", &sr.filt.minMQ, "Minimum mapping quality to consider (lower MQ will be ignored)")
    LONG_INT_PARAM("min-TD", &minTD, "Minimum distance to the tail (lower will be ignored)")
    LONG_INT_PARAM("excl-flag", &sr.filt.exclude_flag, "SAM/BAM FLAGs to be excluded")

    LONG_PARAM_GROUP("Cell/droplet filtering options", NULL)
    LONG_STRING_PARAM("group-list",&groupList, "List of tag readgroup/cell barcode to consider in this run. All other barcodes will be ignored. This is useful for parallelized run")    
    LONG_INT_PARAM("min-total", &minTotalReads, "Minimum number of total reads for a droplet/cell to be considered")
    LONG_INT_PARAM("min-uniq", &minUniqReads, "Minimum number of unique reads (determined by UMI/SNP pair) for a droplet/cell to be considered")
    LONG_INT_PARAM("min-snp", &minCoveredSNPs, "Minimum number of SNPs with coverage for a droplet/cell to be considered")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  std::set<std::string> bcdSet;
  if ( !groupList.empty() ) {
    tsv_reader tsv_group_list(groupList.c_str());
    while( tsv_group_list.read_line() > 0 ) {
      bcdSet.insert(tsv_group_list.str_field_at(0));
    }
    notice("Finished loading %u droplet/cell barcodes to consider", bcdSet.size());
  }

  for(int32_t i=0; i < (int32_t)smIDs.size(); ++i) {
    vr.add_specified_sample(smIDs[i].c_str());
  }

  vr.unlimited_buffer = true;
  vr.vfilt.maxAlleles = 2;
  sr.set_buffer_size(1);
  //sr.unlimited_buffer = true;
  
  vr.init_params();
  sr.init_params();

  int32_t n_warning_no_gtag = 0;
  int32_t n_warning_no_utag = 0;  
  
  if ( outPrefix.empty() )
    error("[E:%s:%d %s] --out parameter is missing",__FILE__,__LINE__,__PRETTY_FUNCTION__);

  char gtag[2] = {0,0};
  char utag[2] = {0,0};    

  if ( tagGroup.empty() ) { // do nothing
  }
  else if ( tagGroup.size() == 2 ) {
    gtag[0] = tagGroup.at(0);
    gtag[1] = tagGroup.at(1);    
  }
  else {
    error("[E:%s:%d %s] Cannot recognize group tag %s. It is suppose to be a length 2 string",__FILE__,__LINE__,__FUNCTION__,tagGroup.c_str());
  }

  if ( tagUMI.empty() ) { // do nothing
  }
  else if ( tagUMI.size() == 2 ) {
    utag[0] = tagUMI.at(0);
    utag[1] = tagUMI.at(1);    
  }
  else {
    error("[E:%s:%d %s] Cannot recognize UMI tag %s. It is suppose to be a length 2 string",__FILE__,__LINE__,__FUNCTION__,tagUMI.c_str());
  }    

  // scan VCF and CRAM simultaneously
  // read a variant first

  sc_dropseq_lib_t scl;

  std::vector<int32_t> snpids;
  //std::vector<int32_t> cellids;  
  
  if ( !vr.read() )
    error("[E:%s Cannot read any single variant from %s]", __PRETTY_FUNCTION__, vr.bcf_file_name.c_str());
  
  if ( !vr.parse_posteriors(vr.cdr.hdr, vr.cursor(), field.c_str(), 0.01) )
    error("[E:%s] Cannot parse posterior probability at %s:%d", __PRETTY_FUNCTION__, bcf_hdr_id2name(vr.cdr.hdr,vr.cursor()->rid), vr.cursor()->pos+1);


  // check if the chromosome names are in the same order between BCF and SAM
  std::map<int32_t,int32_t> rid2tids;
  std::map<int32_t,int32_t> tid2rids;
  std::vector<std::string> tchroms;
  std::vector<std::string> rchroms;
  int32_t ntids = bam_hdr_get_n_targets(sr.hdr);
  int32_t prevrid = -1;
  for(int32_t i=0; i < ntids; ++i) {
    const char* chrom = bam_get_chromi(sr.hdr, i);
    tchroms.push_back(chrom);
    int32_t rid = bcf_hdr_name2id(vr.cdr.hdr, chrom);
    if ( rid >= 0 ) {
      if ( prevrid >= rid ) {
	const char* prevchrom = bcf_hdr_id2name(vr.cdr.hdr, prevrid);
	error("[E:%s] Your VCF/BCF files and SAM/BAM/CRAM files have different ordering of chromosomes. SAM/BAM/CRAM file has %s before %s, but VCF/BCF file has %s after %s", __PRETTY_FUNCTION__, prevchrom, chrom, prevchrom, chrom);
      }
      rid2tids[rid] = i;
      tid2rids[i] = rid;
      rchroms.resize(rid+1);
      rchroms[rid] = chrom;
      prevrid = rid;
    }
  }

  if ( rid2tids.empty() || tid2rids.empty() || ( rid2tids.size() != tid2rids.size() ) ) {
    error("[E:%s] Your VCF/BCF files and SAM/BAM/CRAM files does not have any matching chromosomes, or some chromosome names are duplicated");
  }
  
  int32_t nv = vr.get_nsamples();
  double* gps = new double[nv*3];
  for(int32_t i=0; i < nv * 3; ++i) 
    gps[i] = vr.get_posterior_at(i);
  int32_t snpid = scl.add_snp( vr.cursor()->rid, vr.cursor()->pos+1, vr.cursor()->d.allele[0][0], vr.cursor()->d.allele[1][0], vr.get_af(1), gps);
  snpids.push_back(snpid);

  int32_t ibeg = 0;
  char base, qual;
  int32_t rpos;
  kstring_t readseq = {0,0,0};
  kstring_t readqual = {0,0,0};

  int32_t nReadsMultiSNPs = 0, nReadsSkipBCD = 0, nReadsPass = 0, nReadsRedundant = 0, nReadsN = 0, nReadsLQ = 0, nReadsTMP = 0;
    
  while( sr.read() ) { // read SAM file
    int32_t endpos = bam_endpos(sr.cursor());
    int32_t tid2rid = bcf_hdr_name2id(vr.cdr.hdr, bam_get_chrom(sr.hdr, sr.cursor()));
    if ( tid2rid < 0 ) { // no matching BCF entry in the chromosome, skip;
      continue;
    }

    int32_t n_cleared = vr.clear_buffer_before( bcf_hdr_id2name(vr.cdr.hdr, vr.cursor()->rid), sr.cursor()->core.pos );
    //for(int32_t i=ibeg; i < ibeg+n_cleared; ++i) {
    //  v_umis.clear();
    //}
    ibeg += n_cleared;

    // add new snps
    while( ( !vr.eof ) && ( ( vr.cursor()->rid < tid2rid ) || ( ( vr.cursor()->rid == tid2rid ) && ( vr.cursor()->pos < endpos ) ) ) ) {
      if ( vr.read() ) {
	if ( !vr.parse_posteriors(vr.cdr.hdr, vr.cursor(), field.c_str(), 0.01) )
	  error("[E:%s] Cannot parse posterior probability at %s:%d", __PRETTY_FUNCTION__, bcf_hdr_id2name(vr.cdr.hdr,vr.cursor()->rid), vr.cursor()->pos+1);
	
	gps = new double[nv*3];
	for(int32_t i=0; i < nv * 3; ++i) {
	  gps[i] = vr.get_posterior_at(i);
	}
	snpid = scl.add_snp( vr.cursor()->rid, vr.cursor()->pos+1, vr.cursor()->d.allele[0][0], vr.cursor()->d.allele[1][0], vr.get_af(1), gps);
	snpids.push_back(snpid);
      }
      else {
	//error("Cannot read new SNP");
      }
    }

    // get barcode
    int32_t ibcd = 0;
    if ( tagGroup.empty() ) {
      ibcd = scl.add_cell(".");
    }
    else {
      uint8_t *bcd = (*gtag) ? (uint8_t*) bam_aux_get(sr.cursor(), gtag) : NULL;
      const char* sbcd = ".";
      if ( ( bcd != NULL ) && ( *bcd == 'Z' ) ) {
	sbcd = bam_aux2Z(bcd);
      }
      else {
	if ( n_warning_no_gtag < 10 ) {
	  notice("WARNING: Cannot find Droplet/Cell tag %s from %d-th read %s at %s:%d-%d. Treating all of them as a single group", tagUMI.c_str(), sr.n_read, bam_get_qname(sr.cursor()), bam_get_chrom(sr.hdr, sr.cursor()), sr.cursor()->core.pos, bam_endpos(sr.cursor()));
	}
	else if ( n_warning_no_gtag == 10 ) {
	  notice("WARNING: Suppressing 10+ missing Droplet/Cell tag warnings...");
	}
	++n_warning_no_gtag;	
      }
      
      if ( bcdSet.empty() || ( bcdSet.find(sbcd) != bcdSet.end() ) ) {
	ibcd = scl.add_cell(sbcd);
      }
      else {
	++nReadsSkipBCD;
	continue;
      }
    }

    ++nReadsTMP;

    // get UMI
    std::string sumi(".");
    if ( tagUMI.empty() ) {
      catprintf(sumi,"%x",rand()); // give a random UMI
    }
    else {
      uint8_t *umi = (*utag) ? (uint8_t*) bam_aux_get(sr.cursor(), utag) : NULL;      
      if ( ( umi != NULL ) && ( *umi == 'Z' ) ) {
	sumi = bam_aux2Z(umi);
      }
      else {
	if ( n_warning_no_utag < 10 ) {
	  notice("WARNING: Cannot find UMI tag %s from %d-th read %s at %s:%d-%d. Treating all of them as a single UMI", tagUMI.c_str(), sr.n_read, bam_get_qname(sr.cursor()), bam_get_chrom(sr.hdr, sr.cursor()), sr.cursor()->core.pos, bam_endpos(sr.cursor()));
	}
	else if ( n_warning_no_utag == 10 ) {
	  notice("WARNING: Suppressing 10+ UMI warnings...");
	}
	++n_warning_no_utag;
	//error("[E:%s] Cannot find UMI tag %d %d %x %s %s %x", __PRETTY_FUNCTION__, sr.nbuf, sr.ridx, sr.cursor(), bcd, utag, umi);
      }
    }

    ++scl.cell_totl_reads[ibcd];    
    
    // genotype all reads together
    int32_t nv_pass = 0;
    int32_t nv_redundant = 0;
    int32_t nv_valid = 0;
    int32_t allele, bq;

    //if ( rand() % 10000 == 0 )
    //notice("Reading between %s:%d-%d at %s:%d to %d i=beg=%d, nbuf=%d, vidx=%d, size=%u prevpos=%d", bam_get_chrom(sr.hdr, sr.cursor()), sr.cursor()->core.pos+1, bam_endpos(sr.cursor()), bcf_hdr_id2name(vr.cdr.hdr, scl.snps[ibeg].rid), scl.snps[ibeg].pos, scl.snps[ibeg+vr.nbuf-1].pos, ibeg, vr.nbuf, vr.vidx, vr.vbufs.size(), scl.snps[ibeg-1].pos);
    
    for(int32_t i=ibeg; i < ibeg+vr.nbuf; ++i) {
      bam1_t* b = sr.cursor();
      bam_get_base_and_qual_and_read_and_qual(b, (uint32_t)scl.snps[i].pos-1, base, qual, rpos, &readseq, &readqual);
      if ( rpos == BAM_READ_INDEX_NA ) {
	//if ( rand() % 1000 == 0 ) 
	//notice("Cannot find any informative read between %s:%d-%d at %s:%d", bam_get_chrom(sr.hdr, b), b->core.pos+1, bam_endpos(b), bcf_hdr_id2name(vr.cdr.hdr, scl.snps[i].rid), scl.snps[i].pos);
	continue;
      }
      if ( base == 'N' ) continue;

      ++nv_valid;

      if ( qual-33 < minBQ ) { continue; }
      if ( rpos < minTD-1 ) { continue; }
      if ( rpos + minTD > b->core.l_qseq ) { continue; }

      allele = ( base == scl.snps[i].ref ) ? 0 : ( ( base == scl.snps[i].alt ) ? 1 : 2 );
      bq = qual-33 > capBQ ? capBQ : qual-33;

      if ( scl.add_read(snpids[i], ibcd, sumi.c_str(), allele, bq) )
	++nv_pass;
      else
	++nv_redundant;
    }

    if ( nv_pass > 1 ) ++nReadsMultiSNPs;
    if ( nv_pass > 0 ) ++nReadsPass;
    else if ( nv_redundant > 0 ) ++nReadsRedundant;
    else if ( nv_valid > 0 ) ++nReadsLQ;
    else ++nReadsN;
  }

  if ( n_warning_no_utag > 10 ) 
    notice("WARNING: Suppressed a total of %d UMI warnings...", n_warning_no_utag);
    
  if ( n_warning_no_gtag > 10 ) 
    notice("WARNING: Suppressed a total of %d droplet/cell barcode warnings...", n_warning_no_gtag);
  
  notice("Finished reading %d markers from the VCF file", (int32_t)snpids.size());

  //notice("Finished processing %d reads across %d variants across %d barcodes", nReadsPass, (int32_t)v_poss.size(), (int32_t)bcMap.size(), (int32_t)bcMap.size());
  notice("Total number input reads : %d", sr.n_read);
  notice("Total number of read-QC-passed reads : %d ", sr.n_read - sr.n_skip); //, nReadsN + nReadsUnique + nReadsLQ + nReadsPass);
  notice("Total number of skipped reads with ignored barcodes : %d", nReadsSkipBCD);
  notice("Total number of non-skipped reads with considered barcodes : %d", nReadsTMP);  
  notice("Total number of gapped/noninformative reads : %d", nReadsN);
  notice("Total number of base-QC-failed reads : %d", nReadsLQ);  
  notice("Total number of redundant reads : %d", nReadsRedundant);
  notice("Total number of pass-filtered reads : %d", nReadsPass);
  notice("Total number of pass-filtered reads overlapping with multiple SNPs : %d", nReadsMultiSNPs);

  sr.close();
  //vr.close();

  notice("Starting to prune out cells with too few reads...");
  int32_t nRemoved = 0;
  if ( minTotalReads + minUniqReads + minCoveredSNPs < 0 ) {
    for(int32_t i=0; i < scl.nbcs; ++i) {
      if ( ( scl.cell_totl_reads[i] < minTotalReads ) || ( scl.cell_uniq_reads[i] < minUniqReads) || ( (int32_t)scl.cell_umis[i].size() < minCoveredSNPs ) ) {
	for(std::map<int32_t,sc_snp_droplet_t*>::iterator it = scl.cell_umis[i].begin();
	    it != scl.cell_umis[i].end(); ++it) {
	  delete it->second;
	  scl.snp_umis[it->first].erase(i);
	}
	scl.cell_umis[i].clear();
	++nRemoved;
      }
    }
  }
  notice("Finishing pruning out %d cells with too few reads...", nRemoved);


  if ( (int32_t)snpids.size() != scl.nsnps )
    error("[E:%s snpids.size() = %u != scl.nsnps = %d",__PRETTY_FUNCTION__, snpids.size(), scl.nsnps);

  // Calculate average genotype probability
  double* gp0s = (double*) calloc(scl.nsnps * 3, sizeof(double)); 
  for(int32_t i=0; i < scl.nsnps; ++i) {
    for(int32_t j=0; j < nv; ++j) {
      gp0s[i*3] += scl.snps[i].gps[3*j];
      gp0s[i*3+1] += scl.snps[i].gps[3*j+1];
      gp0s[i*3+2] += scl.snps[i].gps[3*j+2];      
    }
    gp0s[i*3] /= nv;
    gp0s[i*3+1] /= nv;
    gp0s[i*3+2] /= nv;    
  }

  std::map<int64_t, size_t> snpcell2idx;
  std::map<int64_t, size_t> cellsnp2idx;  
  //std::vector<double> gls;

  htsFile* wC = hts_open((outPrefix+".cel.gz").c_str(),"wg");
  htsFile* wV = hts_open((outPrefix+".var.gz").c_str(),"wg");
  htsFile* wP = hts_open((outPrefix+".plp.gz").c_str(),"wg");  
  
  if ( ( wC == NULL ) || ( wV == NULL ) || ( wP == NULL ) )
    error("[E:%s:%d %s] Cannot create %s.* file",__FILE__,__LINE__,__FUNCTION__,outPrefix.c_str());

  notice("Writing cell information");
  std::vector<std::string> v_bcs(scl.bc_map.size());
  for(std::map<std::string,int32_t>::iterator it = scl.bc_map.begin(); it != scl.bc_map.end(); ++it) {
    if ( !v_bcs[it->second].empty() )
      error("Duplicate position for barcode %s at %d", it->first.c_str(), it->second);
    v_bcs[it->second] = it->first;
    //hprintf(wC, "%s\t%d\n", it->first.c_str(), it->second);
  }
  for(int32_t i=0; i < (int32_t)v_bcs.size(); ++i) {
    hprintf(wC, "%d\t%s\t%d\n", i, v_bcs[i].c_str());
  }
  v_bcs.clear();
  notice("Finished writing cell information");
  hts_close(wC);

  notice("Build indices for sparse matrix");
  //size_t k;
  int32_t i;
  //double tmp;
  //for(i=0, k=0; i < scl.nsnps; ++i) {
  for(i=0; i < scl.nsnps; ++i) {    
    hprintf(wV, "%d\t%s\t%d\t%c\t%c\t%.5lf\n", i, rchroms[scl.snps[i].rid].c_str(), scl.snps[i].pos, scl.snps[i].ref, scl.snps[i].alt, 0.5*gp0s[i*3+1] + gp0s[i*3+2]);
    
    std::map<int32_t,sc_snp_droplet_t*>& cells = scl.snp_umis[i];
    if ( cells.empty() ) continue;
    std::map<int32_t,sc_snp_droplet_t*>::iterator it;

    //gls.resize(3*(k+cells.size()));
    //double GLs[3];
    for(it = cells.begin(); it != cells.end(); ++it) {
      //GLs[0] = GLs[1] = GLs[2] = 1.0;
      std::string s_als;
      std::string s_bqs;

      for(sc_snp_droplet_it_t it2=it->second->begin(); it2 != it->second->end(); ++it2) {
	uint8_t al = ( it2->second >> 24 ) & 0x00ff;
	uint8_t bq = ( it2->second >> 16 ) & 0x00ff;
	s_als += ('0' + al);
	s_bqs += (33 + bq);
      }
      hprintf(wP, "%d\t%d\t%s\t%s\n", it->first, i, s_als.c_str(), s_bqs.c_str());
    }
  }
 
  notice("Finished builing indices for sparse matrix");
  hts_close(wP);
  hts_close(wV);  
  
  return 0;
}
