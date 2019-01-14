#include "cramore.h"
#include "sam_filtered_reader.h"
#include "sam_ordered_writer.h"
#include "utils.h"
#include <functional>

typedef std::map<std::string, std::pair<std::string,double> > mux_map_t;
typedef std::map<std::string, std::pair<std::string,double> >::iterator mux_map_it_t;

int32_t cmdCramSimuxlet(int32_t argc, char** argv) {
  SAMFilteredReader sr;
  //sr.verbose = 100000000;
  std::string mux_map_file;
  std::string outPrefix;
  std::string tagGroup = "CB";
  std::string tagUMI = "UB";
  bool drop_other;
  int32_t seed = 0;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&sr.sam_file_name, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("tag-group",&tagGroup, "Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB")
    LONG_STRING_PARAM("tag-UMI",&tagUMI, "Tag representing UMIs. For 10x genomiucs, use UB")

    LONG_PARAM_GROUP("Options for input multiplex file", NULL)
    LONG_STRING_PARAM("mux-map",&mux_map_file,"File containing [BARCODE1] [BARCODE2] [ALPHA1=0.5] [ALPHA2=1-ALPHA1] at each line")
    LONG_INT_PARAM("seed",&seed,"Seed for the UMI hash function")
    LONG_PARAM("drop-other",&drop_other,"Do not create another cell for the remaining pair")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Output file name")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  sr.set_buffer_size(1);
  sr.init_params();

  if ( seed == 0 ) seed = (int32_t)time(NULL);

  if ( mux_map_file.empty() ) error("Missing required option --mux-map");
  if ( outPrefix.empty() ) error("Missing required option --out");

  tsv_reader tsv_mux(mux_map_file.c_str());
  mux_map_t mux_map;
  
  while( tsv_mux.read_line() > 0 ) {
    double alpha = tsv_mux.nfields > 2 ? tsv_mux.double_field_at(2) : 0.5;
    double beta = tsv_mux.nfields > 3 ? tsv_mux.double_field_at(3) : 1-alpha;
    std::string s1(tsv_mux.str_field_at(0));
    std::string s2(tsv_mux.str_field_at(1));

    if ( ( mux_map.find(s1) == mux_map.end() ) && ( mux_map.find(s2) == mux_map.end() ) ) {
      mux_map[s1] = std::make_pair(s2,alpha);
      mux_map[s2] = std::make_pair(s1,beta);
    }
    else {
      error("[E:%s] Barcode %s or %s is appearing multiple times",__PRETTY_FUNCTION__,s1.c_str(),s2.c_str());
    }
  }

  char gtag[2] = {0,0};
  char utag[2] = {0,0};    

  if ( tagGroup.size() == 2 ) {
    gtag[0] = tagGroup.at(0);
    gtag[1] = tagGroup.at(1);    
  }
  else {
    error("[E:%s:%d %s] Cannot recognize group tag %s. It is suppose to be a length 2 string",__FILE__,__LINE__,__FUNCTION__,tagGroup.c_str());
  }

  if ( tagUMI.size() == 2 ) {
    utag[0] = tagUMI.at(0);
    utag[1] = tagUMI.at(1);    
  }
  else {
    error("[E:%s:%d %s] Cannot recognize UMI tag %s. It is suppose to be a length 2 string",__FILE__,__LINE__,__FUNCTION__,tagUMI.c_str());
  }    

  // scan VCF and CRAM simultaneously
  // read a variant first
  SAMOrderedWriter sow(outPrefix.c_str());
  sow.set_hdr(sr.hdr);
  sow.write_hdr();

  std::map<std::string, std::pair<int32_t,int32_t> > bcd_counts;

  while( sr.read() ) { // read SAM file
    bam1_t* b = sr.cursor();
    uint8_t *bcd = (*gtag) ? (uint8_t*) bam_aux_get(b, gtag) : NULL;
    if ( ( bcd != NULL ) && ( *bcd == 'Z' ) ) { // if barcode presents
      std::string sbcd(bam_aux2Z(bcd));
      if ( mux_map.find(sbcd) != mux_map.end() ) { // barcode needs to be multiplexed
	std::pair<std::string,double>& mux_pair = mux_map[sbcd];

	// get UMI
	std::string sumi; // get umi
	uint8_t *umi = (*utag) ? (uint8_t*) bam_aux_get(b, utag) : NULL;
	if ( ( umi != NULL ) && ( *umi == 'Z' ) ) {
	  sumi = bam_aux2Z(umi);
	}
	else {
	  catprintf(sumi, "%x",rand());
	}

	// hash UMI
	uint16_t humi = (str_hash(sumi.c_str()) % UINT16_MAX);
	double dumi = (double)(humi+0.5) / (double)UINT16_MAX;
	std::string new_bcd;

	if ( dumi < mux_pair.second ) { // need to create a new barcode
	  if ( sbcd < mux_pair.first ) {
	    new_bcd = sbcd + "." + mux_pair.first;
	    ++(bcd_counts[new_bcd].first);
	  }
	  else {
	    new_bcd = mux_pair.first + "." + sbcd;
	    ++(bcd_counts[new_bcd].second);	    
	  }
	}
	else if ( !drop_other ) {
	  if ( sbcd > mux_pair.first ) {
	    new_bcd = sbcd + "." + mux_pair.first;
	    ++(bcd_counts[new_bcd].first);
	  }
	  else {
	    new_bcd = mux_pair.first + "." + sbcd;
	    ++(bcd_counts[new_bcd].second);	    
	  }
	}

	if ( !new_bcd.empty() ) { // change the barcode
	  bam_aux_update_str(b, gtag, (int32_t)new_bcd.size()+1, new_bcd.c_str());
	}
      }
    }

    sow.write(b);
  }
  notice("Closing the file");
  
  sr.close();
  sow.close();

  notice("Finished writing BAM file and start writing count file");

  std::map<std::string,std::pair<int32_t,int32_t> >::iterator it;
  htsFile* wf = hts_open((outPrefix + ".mux.counts.txt").c_str(),"w");
  for(it = bcd_counts.begin(); it != bcd_counts.end(); ++it) {
    hprintf(wf, "%s\t%d\t%d\n", it->first.c_str(), it->second.first, it->second.second);
  }
  hts_close(wf);

  notice("Finished writing count files");
  
  return 0;
}
