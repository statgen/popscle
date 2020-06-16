#include "cramore.h"
#include "sc_drop_seq.h"
#include "genomeLoci.h"
#include <sys/stat.h>

// TODO : Record strand info
// TODO : Reduce memory footprint
int32_t cmdCramDscDump2Pileup(int32_t argc, char** argv) {
  std::string inPrefix;
  std::string outPrefix;
  int32_t chunk_id = -1;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Required pptions for input/output files", NULL)
    LONG_STRING_PARAM("in",&inPrefix, "Input prefix from cramore/popscle dsc-dump")
    LONG_STRING_PARAM("out",&outPrefix, "Output prefix compatible with cramore/popscle dsc-pileup")

    LONG_PARAM_GROUP("Additional options for input files", NULL)    
    LONG_INT_PARAM("chunk-id",&chunk_id, "Chunk ID to focus on. By default, it runs across all chunks")    
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( inPrefix.empty() )
    error("Missing required option --in-prefix");

  notice("Reading %s.barcodes.tsv.gz", inPrefix.c_str());
  // read the input barcodes and variants first
  tsv_reader tsv_bcds( (inPrefix + ".barcodes.tsv.gz").c_str() );
  tsv_bcds.read_line(); // skip the first line
  std::vector<std::string> bcds;
  std::vector< std::vector<int> > ibcd2chunks;
  while( tsv_bcds.read_line() > 0 ) {
    if ( tsv_bcds.nfields != 3 )
      error("Cannot parse the line %s because it does not have 3 tokens as expected", tsv_bcds.str.s);
    if ( tsv_bcds.int_field_at(0) != (int32_t)bcds.size() )
      error("Cannot parse the line %s because expected ID is %d but observed is %d", (int32_t)bcds.size(), tsv_bcds.int_field_at(0));
    if ( ibcd2chunks.size() <= tsv_bcds.int_field_at(2) )
      ibcd2chunks.resize(tsv_bcds.int_field_at(2)+1);
    ibcd2chunks[tsv_bcds.int_field_at(2)].push_back(tsv_bcds.int_field_at(0));
    bcds.push_back(tsv_bcds.str_field_at(1));
  }

  tsv_reader tsv_vars;
  std::map<std::string,int32_t> chrom2rid;
  std::vector<std::string> chroms;
  std::vector<int32_t> rids;
  std::vector<int32_t> poss;
  std::vector<std::string> refs;
  std::vector<std::string> alts;
  std::vector<double> afs;
  bool skipVars = false;

  notice("Reading %s.variants.tsv.gz", inPrefix.c_str());  
  if ( tsv_vars.open( (inPrefix + ".variants.tsv.gz").c_str() ) ) {
    htsFile* wV = hts_open((outPrefix+".var.gz").c_str(),"wg");
    hprintf(wV, "#SNP_ID\tCHROM\tPOS\tREF\tALT\tAF\n");
    tsv_vars.read_line();
    while( tsv_vars.read_line() > 0 ) {
      if ( tsv_vars.int_field_at(0) != rids.size() )
	error("Cannot parse the line because expected ID is %d but observed is %d", (int32_t)rids.size(), tsv_vars.int_field_at(0));	
      std::string chr = tsv_vars.str_field_at(1);
      if ( chrom2rid.find(chr) == chrom2rid.end() ) {
	chrom2rid[chr] = chroms.size();
	chroms.push_back(chr);
      }
      rids.push_back(chrom2rid[chr]);
      poss.push_back(tsv_vars.int_field_at(2));
      refs.push_back(tsv_vars.str_field_at(3));
      alts.push_back(tsv_vars.str_field_at(4));
      afs.push_back(tsv_vars.double_field_at(5));

      //error("%d %s %s %.5lf",tsv_vars.int_field_at(2),tsv_vars.str_field_at(3),tsv_vars.str_field_at(4),tsv_vars.double_field_at(5));
      
      hprintf(wV, "%d\t%s\t%d\t%c\t%c\t%.5lf\n", (int32_t)rids.size()-1, chr.c_str(), poss.back(), refs.back()[0], alts.back()[0], afs.back());
    }
    hts_close(wV);
  }
  else {
    notice("No variant file found at %s.variants.tsv.gz. Skipping reading..", inPrefix.c_str());
    skipVars = true;
  }

  // process each chunk separately

  int32_t bcdOffset = 0;
  char buf[65535];

  htsFile* wC = hts_open((outPrefix+".cel.gz").c_str(),"wg");
  htsFile* wP = skipVars ? NULL : hts_open((outPrefix+".plp.gz").c_str(),"wg");
  htsFile* wU = hts_open((outPrefix+".umi.gz").c_str(),"wg");

  hprintf(wC, "#DROPLET_ID\tBARCODE\tNUM.READ\tNUM.UMI\tNUM.UMIwSNP\tNUM.SNP\n");
  if ( !skipVars )
    hprintf(wP, "#DROPLET_ID\tSNP_ID\tALLELES\tBASEQS\n");  
  
  for(int32_t i=0; i < (int32_t)ibcd2chunks.size(); ++i) {
    notice("Processing chunk %d..",i);
    sc_dropseq_lib_t* pscl = new sc_dropseq_lib_t; // for pileup info
    
    for(int32_t j=0; j < (int32_t)rids.size(); ++j) {
      pscl->add_snp(rids[j], poss[j], refs[j][0], alts[j][0], afs[j], NULL);
    }
    std::map<int32_t,int32_t> ibcd_map;
    for(int32_t j=0; j < (int32_t)ibcd2chunks[i].size(); ++j) {
      int32_t ibcd = pscl->add_cell(bcds[ibcd2chunks[i][j]].c_str());
      if ( ibcd != j )
	error("Unexpected barcode ID %d when expecting %d", ibcd, j);
      ibcd_map[ibcd2chunks[i][j]] = ibcd;
    }

    std::vector<std::map<std::string, std::pair<genomeLoci,genomeLoci>> > umiLoci(ibcd2chunks[i].size());  // fwd/rev pair    
    sprintf(buf, "%s.tmp/chunk.%05d.tsv.gz", inPrefix.c_str(), i);
    tsv_reader tsv_chunk(buf);
    while(tsv_chunk.read_line() > 0) {
      int ic = tsv_chunk.int_field_at(0);

      //notice("ic = %d",ic);
      
      std::map<int32_t,int32_t>::iterator it = ibcd_map.find(ic);
      if ( it == ibcd_map.end() )
	error("Cannot find %d from ibcd_map", ic);
      int ib = it->second;
      std::string sumi = tsv_chunk.str_field_at(1);
      bool revStrand = tsv_chunk.str_field_at(2)[0] == '-';
      int nloci = tsv_chunk.int_field_at(3);
      int nvars  = tsv_chunk.int_field_at(4);
      if ( nloci > 0 ) {
	genomeLoci& loci = revStrand ? umiLoci[ib][sumi].second : umiLoci[ib][sumi].first;
	for(int32_t j=0; j < nloci; ++j) {
	  const char* s1 = tsv_chunk.str_field_at(5+j);
	  const char* s2 = strchr(s1,':')+1;
	  const char* s3 = strchr(s2,':')+1;
	  std::string s1_str(s1,s2-s1-1);
	  int32_t s2_int = atoi(s2);
	  int32_t s3_int = atoi(s3);
	  loci.add(s1_str.c_str(), s2_int, s3_int);
	}
	loci.resolveOverlaps();
      }
      if ( nvars > 0 ) {
	for(int32_t j=0; j < nvars; ++j) {
	  const char* s1 = tsv_chunk.str_field_at(5+nloci+j);
	  const char* s2 = strchr(s1,':')+1;
	  const char* s3 = strchr(s2,':')+1;
	  int32_t s1_int = atoi(s1);
	  int32_t s2_int = atoi(s2);
	  int32_t s3_int = atoi(s3);
	  pscl->add_read(s1_int, ib, sumi.c_str(), s2_int, s3_int);
	}
      }
      ++(pscl->cell_totl_reads[ib]);
    }
    
    // write cell barcodes
    for(int32_t j=0; j < (int32_t)ibcd2chunks[i].size(); ++j) {
      hprintf(wC, "%d\t%s\t%d\t%u\t%d\t%u\n", j + bcdOffset, bcds[ibcd2chunks[i][j]].c_str(), pscl->cell_totl_reads[j], umiLoci[j].size(), pscl->cell_uniq_reads[j], pscl->cell_umis[j].size());
    }

    // write UMIs
    for(int32_t j=0; j < (int32_t)umiLoci.size(); ++j) {
      for(std::map<std::string, std::pair<genomeLoci,genomeLoci> >::iterator itu = umiLoci[j].begin();
	  itu != umiLoci[j].end(); ++itu) {
	if ( !( itu->second.first.empty() && itu->second.second.empty() ) ) { // at least fwd or rev is not empty
	  hprintf(wU, "%d\t%s\t%u\t%u", j + bcdOffset, itu->first.c_str(), itu->second.first.totalLength(), itu->second.second.totalLength());
	  for(itu->second.first.rewind(); !itu->second.first.isend(); itu->second.first.next()) {
	    hprintf(wU, "\t%s:%d:%d:+", itu->second.first.it->chrom.c_str(), itu->second.first.it->beg1, itu->second.first.it->end0 - itu->second.first.it->beg1 + 1);
	  }
	  for(itu->second.second.rewind(); !itu->second.second.isend(); itu->second.second.next()) {
	    hprintf(wU, "\t%s:%d:%d:-", itu->second.second.it->chrom.c_str(), itu->second.second.it->beg1, itu->second.second.it->end0 - itu->second.second.it->beg1 + 1);
	  }	  
	  hprintf(wU, "\n");
	}      
      }
    }

    // write pileups
    if ( !skipVars ) {
      for(int32_t j=0; j < pscl->nsnps; ++j) {    
	std::map<int32_t,sc_snp_droplet_t*>& cells = pscl->snp_umis[j];
	if ( cells.empty() ) continue;
	std::map<int32_t,sc_snp_droplet_t*>::iterator it;
	for(it = cells.begin(); it != cells.end(); ++it) {
	  std::string s_als;
	  std::string s_bqs;

	  for(sc_snp_droplet_it_t it2=it->second->begin(); it2 != it->second->end(); ++it2) {
	    uint8_t al = ( it2->second >> 24 ) & 0x00ff;
	    uint8_t bq = ( it2->second >> 16 ) & 0x00ff;
	    s_als += ('0' + al);
	    s_bqs += (33 + bq);
	  }
	  hprintf(wP, "%d\t%d\t%s\t%s\n", it->first + bcdOffset, j, s_als.c_str(), s_bqs.c_str());
	} 
      }
    }
    
    bcdOffset += (int32_t)ibcd2chunks[i].size();

    umiLoci.clear();
    delete pscl;
  }
  if ( !skipVars ) hts_close(wP);
  hts_close(wC);
  hts_close(wU);

  notice("Analysis finished");

  return 0;
}
