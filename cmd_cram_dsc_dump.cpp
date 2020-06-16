#include "cramore.h"
#include "bcf_filtered_reader.h"
#include "sam_filtered_reader.h"
#include "sc_drop_seq.h"
#include "genomeLoci.h"
#include <sys/stat.h>
#include <cassert>

// TODO : Record strand info
// TODO : Reduce memory footprint
int32_t cmdCramDscDump(int32_t argc, char** argv) {
  // input/output files
  SAMFilteredReader sr;
  BCFFilteredReader vr;
  std::string region;         // process only specific regions
  int32_t nchunks = 100;    // number of cell barcodes to store in a chunk
  uint32_t seed = 0x811c9dc5;
  std::string outPrefix;      // outPrefix.map, outPrefix.tmp/ will be created
  
  // tags for droplet barcodes and UMI barcodes
  std::string tagGroup("CB");
  std::string tagUMI("UB");
  int32_t capBQ = 40;
  int32_t minBQ = 13;
  int32_t minTD = 0;
  sr.filt.exclude_flag = 0x0f04; // read-level filter to exclude reads
  sr.filt.minMQ = 20;            // filter to cap the mapping quality 
  
  // input options for VCFs (sites)
  std::vector<std::string> smIDs;
  vr.vfilt.maxAlleles = 2;
  
  // output options
  vr.verbose = 10000;
  sr.verbose = 1000000;
  // options to filter droplets
  std::string groupList;     // (optional) focuses on specific barcodes 
  bool skipUmiFlag = false;  //  skip writing UMI info. Only stores variant-focused pileups
  bool skipEmptyGroup = false; // skip reads with no cell barcodes without assiging '.' as barcode
  bool skipEmptyUMI = false;   // skip reads with no UMI information without assinging random UMI

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&sr.sam_file_name, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("region",&region, "Region to be focused on")    
    LONG_STRING_PARAM("tag-group",&tagGroup, "Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB")
    LONG_STRING_PARAM("tag-UMI",&tagUMI, "Tag representing UMIs. For 10x genomiucs, use UB")
    LONG_INT_PARAM("excl-flag", &sr.filt.exclude_flag, "SAM/BAM FLAGs to be excluded")    

    LONG_PARAM_GROUP("Options for input VCF/BCF", NULL)
    LONG_STRING_PARAM("vcf",&vr.bcf_file_name, "Input VCF/BCF file, containing the AC and AN field")
    LONG_MULTI_STRING_PARAM("sm",&smIDs, "List of sample IDs to compare to (default: use all)")
    LONG_STRING_PARAM("sm-list",&vr.sample_id_list, "File containing the list of sample IDs to compare")        

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Output file prefix")
    LONG_INT_PARAM("sam-verbose",&sr.verbose, "Verbose message frequency for SAM/BAM/CRAM")
    LONG_INT_PARAM("vcf-verbose",&vr.verbose, "Verbose message frequency for VCF/BCF")
    LONG_PARAM("skip-umi", &skipUmiFlag, "Do not generate [prefix].umi.gz file, which stores the regions covered by each barcode/UMI pair")
    LONG_PARAM("skip-empty-group", &skipEmptyGroup, "Skip read that does not have group (e.g. cell barcode) information. By default it assigns barcode '.'")
    LONG_PARAM("skip-empty-umi", &skipEmptyUMI, "Skip read that does not have UMI (e.g. cell barcode) information. By default it assigns all reads as a single UMI. To consider them all independent reads, you need to set --tag-UMI '' (empty)")
    LONG_INT_PARAM("chunks", &nchunks, "Number of chunks to store barcodes randomly")
    LONG_INT_PARAM("seed", &seed, "Seed for random number generator")

    LONG_PARAM_GROUP("SNP-overlapping Read filtering Options", NULL)
    LONG_INT_PARAM("cap-BQ", &capBQ, "Maximum base quality (higher BQ will be capped)")
    LONG_INT_PARAM("min-BQ", &minBQ, "Minimum base quality to consider (lower BQ will be skipped)")
    LONG_INT_PARAM("min-MQ", &sr.filt.minMQ, "Minimum mapping quality to consider (lower MQ will be ignored)")
    LONG_INT_PARAM("min-TD", &minTD, "Minimum distance to the tail (lower will be ignored)")

    LONG_PARAM_GROUP("Cell/droplet filtering options", NULL)
    LONG_STRING_PARAM("group-list",&groupList, "List of tag readgroup/cell barcode to consider in this run. All other barcodes will be ignored. This is useful for parallelized run")    
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

  if ( !region.empty() ) {
    sr.target_region = region;
    vr.target_region = region;
  }

  bool bcfEmpty = vr.bcf_file_name.empty();
  if ( bcfEmpty ) {
      notice("No VCF input file was specified. Skipping reading..");
  }
  else {  
    notice("Initializing BCF reader..");
    vr.init_params();
  }
  notice("Initializing SAM reader..");  
  sr.init_params();

  int32_t n_warning_no_gtag = 0;
  int32_t n_warning_no_utag = 0;

  std::vector<std::string> bcds;
  std::vector<int32_t>     ichunks;
  std::map<std::string,int32_t> bcd2idx;
  std::vector<int32_t>     var_rids;
  std::vector<int32_t>     var_poss;
  std::vector<char>        var_refs;
  std::vector<char>        var_alts;
  std::vector<double>      var_afs;
  std::map<std::string,int32_t> var2idx;  

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

  //sc_dropseq_lib_t scl;

  char buf[65535];

  // create output directory
  std::string outDir = outPrefix + ".tmp";
  std::string outMap = outPrefix + ".map.gz";
  int32_t ret = mkdir(outDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if ( ret < 0 ) {
    if ( errno == EEXIST ) {
      warning("Warning: The directory %s already exists. Some existing files in the directory may be replaced", outPrefix.c_str());
    }
    else {
      error("Cannot create the directory %s", outPrefix.c_str());
    }
  }
  
  // create chunk files
  std::vector<htsFile*> wfs(nchunks);
  if ( nchunks > 99999 ) error("Cannot use %d (more than 99999) chunks", nchunks);
  else if ( nchunks < 1 ) error("Cannot use %d (less than 1) chunk", nchunks);

  // open the chunk filess simultaneously
  for(int32_t i=0; i < nchunks; ++i) {
    sprintf(buf, "%s/chunk.%05d.tsv.gz", outDir.c_str(), i);
    wfs[i] = hts_open(buf, "wg");
  }

  // check if the chromosome names are in the same order between BCF and SAM, only if region is not empty
  std::map<int32_t,int32_t> rid2tids;
  std::map<int32_t,int32_t> tid2rids;
  std::vector<std::string> tchroms;
  std::vector<std::string> rchroms;
  int32_t ntids = bam_hdr_get_n_targets(sr.hdr);
  int32_t prevrid = -1;
  for(int32_t i=0; i < ntids; ++i) {
    const char* chrom = bam_get_chromi(sr.hdr, i);
    tchroms.push_back(chrom);
    if ( !bcfEmpty ) {
      int32_t rid = bcf_hdr_name2id(vr.cdr.hdr, chrom);
      if ( rid >= 0 ) {
	if ( prevrid >= rid ) {
	  const char* prevchrom = bcf_hdr_id2name(vr.cdr.hdr, prevrid);
	  if ( region.empty() )
	    error("[E:%s] Your VCF/BCF files and SAM/BAM/CRAM files have different ordering of chromosomes. SAM/BAM/CRAM file has %s before %s, but VCF/BCF file has %s after %s", __PRETTY_FUNCTION__, prevchrom, chrom, prevchrom, chrom);
	}
	rid2tids[rid] = i;
	tid2rids[i] = rid;
	rchroms.resize(rid+1);
	rchroms[rid] = chrom;
	prevrid = rid;
      }
    }
  }


  if ( !bcfEmpty ) {
    if ( !vr.read() ) // Should we throw error, or warning?
      error("[E:%s Cannot read any single variant from %s]", __PRETTY_FUNCTION__, vr.bcf_file_name.c_str());
  
    double af = vr.calculate_af(true);
    bcf1_t* v = vr.cursor();
    sprintf(buf, "%d:%d:%c:%c", (int32_t)v->rid, (int32_t)v->pos+1, (char)v->d.allele[0][0], (char)v->d.allele[1][0]);
    std::map<std::string,int32_t>::iterator it = var2idx.find(buf);
    if ( it == var2idx.end() ) {
      var_rids.push_back(v->rid);
      var_poss.push_back(v->pos+1);
      var_refs.push_back(v->d.allele[0][0]);
      var_alts.push_back(v->d.allele[1][0]);
      var_afs.push_back(af);
      var2idx[buf] = (int32_t)var_rids.size()-1;
    }
  }
    
  //if ( rid2tids.empty() || tid2rids.empty() || ( rid2tids.size() != tid2rids.size() ) ) {
  if ( !bcfEmpty && ( rid2tids.empty() || tid2rids.empty() ) ) {
    error("Your VCF/BCF files and SAM/BAM/CRAM files does not have any matching chromosomes, or some chromosome names are duplicated. This usually can be resolved by running 'bcftools reheader -f [ref.fasta.fai] [in.vcf.gz] > [out.vcf.gz]'");
    //error("[E:%s] Your VCF/BCF files and SAM/BAM/CRAM files does not have any matching chromosomes, or some chromosome names are duplicated");
  }

  int32_t ibeg = 0;
  char base, qual;
  int32_t rpos;
  kstring_t readseq = {0,0,0};
  kstring_t readqual = {0,0,0};

  int32_t nReadsMultiSNPs = 0, nReadsSkipBCD = 0, nReadsPass = 0, nReadsRedundant = 0, nReadsN = 0, nReadsLQ = 0, nReadsTMP = 0, nReadsFilt = 0;

  //notice("foo");

  //std::vector<std::map<std::string, std::pair<genomeLoci,genomeLoci> > > umiLoci;  // fwd/rev pair

  while( sr.read() ) { // read SAM file, processing each individual read
    bam1_t* b = sr.cursor();    
    int32_t endpos = bam_endpos(b);
    const char* chrom = bam_get_chrom(sr.hdr, b);
    int32_t tid2rid = bcfEmpty ? -1 : bcf_hdr_name2id(vr.cdr.hdr, chrom);
    bool noBCF = bcfEmpty || vr.eof;

    int32_t bamFlag = bam_get_flag(b);
    if ( bamFlag & sr.filt.exclude_flag )  {
      ++nReadsFilt;
      continue;
    }
    bool revStrand = (bamFlag & 0x0010) ? true : false; // get strand information
    
    if ( tid2rid < 0 ) { // no matching BCF entry in the chromosome, skip;
      noBCF = true;
    }

    if ( ( !noBCF ) && ( vr.cursor()->rid >= rchroms.size() ) ) { // new contig absent in the header was introduced
      error("The contig %s was absent in the VCF header. This happens when contigs in VCF header do not match to actual records. Run 'bcftools reheader -f [ref.fasta.fai] [in.vcf.gz] > [out.vcf.gz] to resolve this problem", chrom);
    }    
    
    int32_t n_cleared = ( noBCF || vr.eof ) ? 0 : vr.clear_buffer_before( bcf_hdr_id2name(vr.cdr.hdr, vr.cursor()->rid), b->core.pos );
    ibeg += n_cleared;
    
    // add new snps
    if ( !noBCF ) {
      while( ( !vr.eof ) && ( ( vr.cursor()->rid < tid2rid ) || ( ( vr.cursor()->rid == tid2rid ) && ( vr.cursor()->pos < endpos ) ) ) ) {
	if ( vr.read() ) {
	  double af = vr.calculate_af(true);
	  bcf1_t* v = vr.cursor();
	  sprintf(buf, "%d:%d:%c:%c", (int32_t)v->rid, (int32_t)v->pos+1, (char)v->d.allele[0][0], (char)v->d.allele[1][0]);
	  std::map<std::string,int32_t>::iterator it = var2idx.find(buf);
	  if ( it == var2idx.end() ) {
	    var_rids.push_back(v->rid);
	    var_poss.push_back(v->pos+1);
	    var_refs.push_back(v->d.allele[0][0]);
	    var_alts.push_back(v->d.allele[1][0]);
	    var_afs.push_back(af);
	    var2idx[buf] = (int32_t)var_rids.size()-1;
	  }
	  //else 
	  //  snpid = it->second;
	}
	else {
	  //error("Cannot read new SNP");
	}
      }
    }
    
    // get barcode
    int32_t ibcd = 0;
    uint8_t *bcd = (*gtag) ? (uint8_t*) bam_aux_get(b, gtag) : NULL;
    const char* sbcd = ".";
    if ( ( bcd != NULL ) && ( *bcd == 'Z' ) ) {
      sbcd = bam_aux2Z(bcd);
    }
    else {
      if ( skipEmptyGroup ) { // skip reads with no cell barcodes
	++nReadsSkipBCD;
	continue;
      }      
      if ( n_warning_no_gtag < 10 ) {
	notice("WARNING: Cannot find Droplet/Cell tag %s from %d-th read %s at %s:%d-%d. Treating all of them as a single group", tagUMI.c_str(), sr.n_read, bam_get_qname(b), chrom, b->core.pos, endpos);
      }
      else if ( n_warning_no_gtag == 10 ) {
	notice("WARNING: Suppressing 10+ missing Droplet/Cell tag warnings...");
      }
      ++n_warning_no_gtag;
    }
    
    // if barcode is a valid barcode
    if ( bcdSet.empty() || ( bcdSet.find(sbcd) != bcdSet.end() ) ) {
      std::map<std::string,int32_t>::iterator it = bcd2idx.find(sbcd);
      if ( it == bcd2idx.end() ) {
	bcds.push_back(sbcd);
	int32_t ichunk = (int32_t)(str_hash(sbcd, seed) % nchunks); // chunk index
	ichunks.push_back(ichunk);
	ibcd = bcd2idx[sbcd] = (int32_t)bcds.size()-1;	  
      }
      else {
	ibcd = it->second;
      }
    }
    else {
      ++nReadsSkipBCD;
      continue;
    }

    //if ( bcd ) free(bcd);

    ++nReadsTMP;

    // get UMI
    std::string sumi(".");
    if ( tagUMI.empty() ) { // if UMI tag was not defined
      sumi.clear();
      catprintf(sumi,"%x",rand()); // give a random UMI
    }
    else {
      uint8_t *umi = (*utag) ? (uint8_t*) bam_aux_get(b, utag) : NULL;      
      if ( ( umi != NULL ) && ( *umi == 'Z' ) ) {
	sumi = bam_aux2Z(umi);
      }
      else {
	if ( skipEmptyUMI ) {
	  ++nReadsSkipBCD;
	  continue;
	}
	if ( n_warning_no_utag < 10 ) {
	  notice("WARNING: Cannot find UMI tag %s from %d-th read %s at %s:%d-%d. Treating all of them as a single UMI", tagUMI.c_str(), sr.n_read, bam_get_qname(b), bam_get_chrom(sr.hdr, b), b->core.pos, endpos);
	}
	else if ( n_warning_no_utag == 10 ) {
	  notice("WARNING: Suppressing 10+ UMI warnings...");
	}
	++n_warning_no_utag;
	//error("[E:%s] Cannot find UMI tag %d %d %x %s %s %x", __PRETTY_FUNCTION__, sr.nbuf, sr.ridx, sr.cursor(), bcd, utag, umi);
      }
      //if ( umi ) free(umi);
    }

    // we just need to print CHROM POS CIGAR STRAND
    genomeLoci loci;
    if ( !skipUmiFlag && b->core.n_cigar ) { // go over CIGAR string to find 'M's
      //int32_t rlen = b->core.l_qseq;
      int32_t cpos = b->core.pos;
      int32_t rpos = 0;
      kstring_t str = {0, 0, 0};
      
      uint32_t *cigar = bam_get_cigar(b);
      for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
	char op = bam_cigar_opchr(cigar[i]);
	str.l = 0;
	kputw(bam_cigar_oplen(cigar[i]), &str);
	char* stop;
	int32_t len = strtol(str.s, &stop, 10);
	assert(stop);
	
	// if M is observed, add the region to loci
	if (op=='M') {
	  loci.add(chrom, cpos+1, cpos+len);
	  cpos += len;
	  rpos += len;
	}
	else if ( ( op=='D' ) || ( op=='N' ) ) {
	  cpos += len;
	}
	else if (op=='S' || op=='I') {
	  rpos += len;
	}
      }
      loci.resolveOverlaps();
      //if ( cigar ) free(cigar);
    }
    
    //if ( ibcd % 1000 == 0 )
    //  notice("%d %s %d", ibcd, bcds[ibcd].c_str(), ichunks[ibcd]);

    // genotype all reads together
    //int32_t nv_pass = 0;
    //int32_t nv_redundant = 0;
    int32_t nv_valid = 0;
    int32_t allele, bq;

    //for(int32_t i=ibeg; i < ibeg+vr.nbuf; ++i) {
    std::vector<std::string> varBQs;
    if ( !noBCF ) {
      for(int32_t i=ibeg; i < (int32_t)var_rids.size(); ++i) {      
	//if ( i >= (int32_t) ) var_rids.size() continue;
	bam_get_base_and_qual_and_read_and_qual(b, (uint32_t)var_poss[i]-1, base, qual, rpos, &readseq, &readqual);
	if ( rpos == BAM_READ_INDEX_NA ) {
	  continue;
	}
	if ( base == 'N' ) continue;
	
	++nv_valid;
	
	if ( qual-33 < minBQ ) { continue; }
	if ( rpos < minTD-1 ) { continue; }
	if ( rpos + minTD > b->core.l_qseq ) { continue; }
	
	allele = ( base == var_refs[i] ) ? 0 : ( ( base == var_alts[i] ) ? 1 : 2 );
	bq = qual-33 > capBQ ? capBQ : qual-33;

	sprintf(buf,"%d:%d:%d", i, allele, bq);
	varBQs.push_back(buf);
	//hprintf(wfs[ichunks[ibcd]],"\t%d:%d:%d", i, allele, bq);
      }
    }

    // write [BARCODE_ID] [UMI] [STRAND] [N.LOCI] [N.VAR] [LOCUS1] .. [LOCUSN] [N.VAR] [VAR1] .. [VARN]
    if ( !skipUmiFlag || !varBQs.empty() ) { // write marker info
      hprintf(wfs[ichunks[ibcd]],"%d\t%s\t%c\t%d\t%d", ibcd, sumi.c_str(), revStrand ? '-' : '+', (int32_t)loci.loci.size(), (int32_t)varBQs.size());
      for(loci.rewind(); !loci.isend(); loci.next()) {
	hprintf(wfs[ichunks[ibcd]],"\t%s:%d:%d", loci.it->chrom.c_str(), loci.it->beg1, loci.it->end0 - loci.it->beg1 + 1);
      }
      for(int32_t i=0; i < (int32_t)varBQs.size(); ++i)
	hprintf(wfs[ichunks[ibcd]],"\t%s", varBQs[i].c_str());
      hprintf(wfs[ichunks[ibcd]],"\n");      
    }
  }

  for(int32_t i=0; i < nchunks; ++i) 
    hts_close(wfs[i]);

  // write barcode and SNP information
  htsFile* wC = hts_open((outPrefix + ".barcodes.tsv.gz").c_str(), "wg");
  hprintf(wC, "ID\tBARCODE\tCHUNK\n" );
  for(int32_t i=0; i < (int32_t)bcds.size(); ++i) {
    hprintf(wC,"%d\t%s\t%d\n", i, bcds[i].c_str(), ichunks[i]);
  }
  hts_close(wC);

  if ( !bcfEmpty ) {
    htsFile* wV = hts_open((outPrefix + ".variants.tsv.gz").c_str(), "wg");
    hprintf(wV, "ID\tCHROM\tPOS\tREF\tALT\tAF\n");
    for(int32_t i=0; i < (int32_t)var_rids.size(); ++i) {
      hprintf(wC,"%d\t%s\t%d\t%c\t%c\t%.5lf\n", i, rchroms[var_rids[i]].c_str(), var_poss[i], var_refs[i], var_alts[i], var_afs[i]);
    }
    hts_close(wV);
  }

  if ( n_warning_no_utag > 10 ) 
    notice("WARNING: Suppressed a total of %d UMI warnings...", n_warning_no_utag);
    
  if ( n_warning_no_gtag > 10 ) 
    notice("WARNING: Suppressed a total of %d droplet/cell barcode warnings...", n_warning_no_gtag);

  /*
  notice("Finished reading %d markers from the VCF file", (int32_t)snpids.size());

  //notice("Finished processing %d reads across %d variants across %d barcodes", nReadsPass, (int32_t)v_poss.size(), (int32_t)bcMap.size(), (int32_t)bcMap.size());
  notice("Total number input reads : %d", sr.n_read);
  notice("Total number of unmapped or QC-failed reads : %d", nReadsFilt);
  notice("Total number of read-QC-passed reads : %d ", sr.n_read - sr.n_skip); //, nReadsN + nReadsUnique + nReadsLQ + nReadsPass);
  notice("Total number of skipped reads with ignored barcodes : %d", nReadsSkipBCD);
  notice("Total number of non-skipped reads with considered barcodes : %d", nReadsTMP);  
  notice("Total number of gapped/noninformative reads : %d", nReadsN);
  notice("Total number of base-QC-failed reads : %d", nReadsLQ);  
  notice("Total number of redundant reads : %d", nReadsRedundant);
  notice("Total number of pass-filtered reads : %d", nReadsPass);
  notice("Total number of pass-filtered reads overlapping with multiple SNPs : %d", nReadsMultiSNPs);
  */

  //sr.close();
  //if ( !bcfEmpty ) vr.close();
  notice("Analysis finished");

  return 0;
}
