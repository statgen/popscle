#include "cramore.h"

#include "dropseq.h"

struct drop_read {
  int32_t tid;
  std::string tname;
  std::string umi;
  std::string barcode;
  std::string seq;
  std::string qual;
  int32_t nh;
  int32_t nm;

  drop_read() : tid(0), nh(0), nm(0) {}
};

typedef struct drop_read drop_read_t;

struct drop_count {
  int32_t umis;
  int32_t reads;
  double fumis;
  int32_t etc;
  double etc2;

  drop_count() : umis(0), reads(0), fumis(0), etc(0), etc2(0) {}
};

typedef struct drop_count drop_count_t;


drop_read_t* bam_drop_read(bam_hdr_t* h, bam1_t* b, const char* NH, const char* NM) {
  uint8_t* s = bam_aux_get(b, NH);
  if ( !s ) {
    error("[E:%s:%d %s] Cannot find %c%c tag in record\n",__FILE__,__LINE__,__FUNCTION__, NH[0], NH[1]);
    return NULL;
  }
  int32_t vNH = bam_aux2i(s);
    
  s = bam_aux_get(b, NM);
  if ( !s ) {
    error("[E:%s:%d %s] Cannot find %c%c tag in record\n",__FILE__,__LINE__,__FUNCTION__, NM[0], NM[1]);
    return NULL;
  }
  int32_t vNM = bam_aux2i(s);
  
  drop_read_t* read = new drop_read_t;
  read->nh = vNH;
  read->nm = vNM;
  
  // extract barcode
  char *prn = bam_get_qname(b);
  char *pbc = NULL;
  char *pumi = prn;
  char *ptmp = NULL;
  while( ( ptmp = strchr(pumi, ':') ) != NULL ) {
    pbc = pumi+1;
    pumi = ptmp+1;
  }

  read->barcode.assign(pbc,pumi-pbc-1);
  std::transform(read->barcode.begin(), read->barcode.end(), read->barcode.begin(), ::toupper);
  
  read->umi.assign(pumi);
  read->tid = b->core.tid;
  read->tname = h->target_name[read->tid];

  return read;
}

int32_t cmdScKallistoCount(int32_t argc, char** argv) {
  std::string inFile;
  std::string outPrefix;  
  std::string tagNH("NH");
  std::string tagNM("NM");  
  int32_t maxNM = 3;
  
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("SAM/BAM/CRAM Input Options", NULL)
    LONG_STRING_PARAM("sam",&inFile, "Input SAM/BAM/CRAM file. Sorted by readnames")
    
    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_STRING_PARAM("tag-NH",&tagNH, "Tag indicating the number of multiple mapping")
    LONG_STRING_PARAM("tag-NM",&tagNM, "Tag indicating the number of mismatches")
    LONG_INT_PARAM("max-NM",&maxNM, "Maximum number of mismatches allowed")    
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outPrefix, "Output prefix")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inFile.empty() || outPrefix.empty()  ) {
    error("[E:%s:%d %s] --sam, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  samFile* in = NULL;
  bam_hdr_t *header = NULL;
  
  if ( ( in = sam_open(inFile.c_str(), "r") ) == 0 ) {
    error("[E:%s:%d %s] Cannot open file %s\n",__FILE__,__LINE__,__FUNCTION__,inFile.c_str());    
  }

  if ( ( header = sam_hdr_read(in) ) == 0 ) {
    error("[E:%s:%d %s] Cannot open header from %s\n",__FILE__,__LINE__,__FUNCTION__,inFile.c_str());
  }

  bam1_t *b = bam_init1();
  int32_t r;

  // try to read all files
  char NH[2]; NH[0] = tagNH[0]; NH[1] = tagNH[1];
  char NM[2]; NM[0] = tagNM[0]; NM[1] = tagNM[1];

  std::vector<drop_read_t*> v_reads;
  std::map<int32_t, std::string> tid2name;

  DropLibrary dl;

  while( ( r = sam_read1(in, header, b) ) >= 0 ) {
    drop_read_t* rd = bam_drop_read(header, b, NH, NM);
    if ( rd == NULL ) error("[E:%s:%d %s] Cannot extract a read from a BAM",__FILE__,__LINE__,__FUNCTION__);
    v_reads.push_back(rd);
    int32_t minNM = rd->nm;
    if ( rd->nh > 1 ) {
      for(int32_t i=1; i < rd->nh; ++i) {
	if ( ( r = sam_read1(in, header, b) ) < 0 ) {
	  error("[E:%s:%d %s] No more read observed despite NH tag",__FILE__,__LINE__,__FUNCTION__);
	}
	drop_read_t* rd2 = bam_drop_read(header, b, NH, NM);
	if ( ( rd->barcode != rd2->barcode ) || ( rd->umi != rd2->umi ) ) {
	  error("[E:%s:%d %s] Barcode or UMI does not match\n",__FILE__,__LINE__,__FUNCTION__);
	}
	if ( minNM > rd2->nm )
	  minNM = rd2->nm;
	v_reads.push_back(rd2);
      }
    }

    // count only the reads with minimal matches
    if ( minNM <= maxNM ) {
      for(size_t i=0; i < v_reads.size(); ++i) {
	rd = v_reads[i];
	if ( rd->nm == minNM ) {
	  //printf("%s\t%s\t%d\t%s\t%d\t%d\n", rd->barcode.c_str(), rd->umi.c_str(), rd->tid, rd->tname.c_str(), rd->nh, rd->nm);
	  dl.addRead( rd->barcode, rd->tid, rd->umi, rd->nh );

	  if ( tid2name.find(rd->tid) == tid2name.end() ) {
	    tid2name[rd->tid] = rd->tname;
	  }
	}
	delete rd;
      }
    }
    else {
      for(size_t i=0; i < v_reads.size(); ++i) {
	delete v_reads[i];
      }
    }
    v_reads.clear();
  }

  FILE* fpCT = fopen((outPrefix+".cell.tx.cnts").c_str(), "w");
  if ( fpCT == NULL )
    error("[E:%s:%d %s] Cannot write file",__FILE__,__LINE__,__FUNCTION__);
  fprintf(fpCT,"#Cell\tTx\tUMIs\tReads\tfUMIs\n");
  
  FILE* fpC = fopen((outPrefix+".cell.cnts").c_str(), "w");
  if ( fpC == NULL )
    error("[E:%s:%d %s] Cannot write file",__FILE__,__LINE__,__FUNCTION__);    
  fprintf(fpC,"#Cell\tUMIs\tReads\tfUMIs\tnTx\tumiENST\t%%ENST\n");

  FILE* fpT = fopen((outPrefix+".tx.cnts").c_str(), "w");
  if ( fpT == NULL )
    error("[E:%s:%d %s] Cannot write file",__FILE__,__LINE__,__FUNCTION__);    
  fprintf(fpT,"#Tx\tUMIs\tReads\tfUMIs\tnCell\n");

  std::map<int32_t,drop_count_t> t2cnt;

  notice("Writing digital expression matrices");

  for(sc_map_it_t itC = dl.mapCell.begin(); itC != dl.mapCell.end(); ++itC) {
    drop_count_t cCnt;
    
    it_map_t& mT = itC->second->mapTranscript;
    for(it_map_it_t itT = mT.begin(); itT != mT.end(); ++itT) {
      sr_map_t& mR = itT->second->mapRead;
      drop_count_t ctCnt;
      for(sr_map_it_t itR = mR.begin(); itR != mR.end(); ++itR) {
	++ctCnt.umis;
	ctCnt.reads += itR->second.first;
	ctCnt.fumis += (1.0/itR->second.second);
      }

      std::string& tname = tid2name[itT->first];

      cCnt.umis += ctCnt.umis;
      cCnt.reads += ctCnt.reads;
      cCnt.fumis += ctCnt.fumis;
      ++cCnt.etc;
      if ( tname.compare(0, 4, "ENST") == 0 ) {
	cCnt.etc2 += ctCnt.fumis;
      }

      drop_count_t& tCnt = t2cnt[itT->first];

      tCnt.umis += ctCnt.umis;
      tCnt.reads += ctCnt.reads;
      tCnt.fumis += ctCnt.fumis;
      ++tCnt.etc;      
      
      fprintf(fpCT, "%s\t%s\t%d\t%d\t%.3lf\n",itC->first.c_str(), tname.c_str(), ctCnt.umis, ctCnt.reads, ctCnt.fumis);
    }
    
    fprintf(fpC, "%s\t\t%d\t%d\t%.3lf\t%d\t%.3lf\t%.3lf\n",itC->first.c_str(), cCnt.umis, cCnt.reads, cCnt.fumis, cCnt.etc, cCnt.etc2, cCnt.etc2/cCnt.fumis);   
  }

  for(std::map<int32_t,drop_count_t>::iterator it = t2cnt.begin(); it != t2cnt.end(); ++it) {
    drop_count_t& tCnt = it->second;    
    fprintf(fpT, "%s\t\t%d\t%d\t%.3lf\t%d\n", tid2name[it->first].c_str(), tCnt.umis, tCnt.reads, tCnt.fumis, tCnt.etc);              
  }

  bam_destroy1(b);
  
  fclose(fpCT);
  fclose(fpC);
  fclose(fpT);  

  return 0;
}

