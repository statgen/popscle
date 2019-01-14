#include "gtf.h"
#include "tsv_reader.h"

gtfGene::gtfGene(const char* _seqname, int32_t _start, int32_t _end,
		 bool _fwdStrand, std::string& _gid, std::string& _gname,
		 std::string& _gtype)
  : gtfElement(_start, _end, "gene", NULL),
    geneId(_gid), geneName(_gname),
    geneType(_gtype), seqname(_seqname),
    fwdStrand(_fwdStrand) {}

gtfTranscript::gtfTranscript(int32_t _start, int32_t _end, std::string& _tid, std::string& _ttype, gtfElement* _parent)
  : gtfElement(_start, _end, "transcript", _parent),
    transcriptId(_tid), transcriptType(_ttype) {}

gtfElement::gtfElement(int32_t _start, int32_t _end, const char* _type, gtfElement* _parent)
  : parent(_parent), type(_type), locus(_start, _end) {}

gtfCDS::gtfCDS(int32_t _start, int32_t _end, const char* sframe, gtfElement* _parent)
  : gtfElement(_start, _end, "CDS", _parent) {
    // check whether frame is valid integer
    if ( sframe[0] == '0' )      { frame = 0; }
    else if ( sframe[0] == '1' ) { frame = 1; }
    else if ( sframe[0] == '2' ) { frame = 2; }
    else
     error("[E:%s:%d:%s] Unrecognized frame string %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, sframe);
}

gtf::gtf(const char* filename, bool proteinCodingOnly, bool addChrPrefix, bool removeChrPrefix) :
  maxGeneLength(0), maxTranscriptLength(0), maxExonLength(0),
  maxCDSLength(0), maxUTRLength(0), maxStartCodonLength(0),
  maxStopCodonLength(0)
{
  tsv_reader tr;
  tr.delimiter = '\t';
  tr.open(filename);
  int nfields = 0;
  // read each line
  notice("Started reading GTF file %s...", filename);
  int32_t line = 0;
  char seqname[65535];
  //int32_t nwarnings = 0;
  for( line = 0; ( nfields = tr.read_line() ) > 0; ++line) {
    if ( line % 1000000 == 0 )
      notice("Reading line %d from %s...", line, filename);

    // Process chromosome name
    const char* tmp_seqname = tr.str_field_at(0);    
    if ( tmp_seqname[0] == '#' ) continue; // ignore meta lines
    if ( addChrPrefix ) {
      if ( strncmp(tmp_seqname,"chr",3) == 0 )
	error("--gtf-add-chr option is on, but chr is already included in the chromosome name");
      seqname[0] = 'c'; seqname[1] = 'h'; seqname[2] = 'r';
      strcpy(seqname+3, tmp_seqname);
    }
    else if ( removeChrPrefix ) {
      if ( strncmp(tmp_seqname,"chr",3) == 0 ) {
	strcpy(seqname, tmp_seqname+3);
      }
      else {
	//if ( nwarnings < 10 ) 
	//warning("--gtf-remove-chr option is on, but the chromosome name does %s not start with chr. Including the entry without chr prefix", tmp_seqname);
	//++nwarnings;	
	//if ( nwarnings == 10 )
	//warning("10 warnings in same kind were observed. Supressing further warnings...");
	strcpy(seqname, tmp_seqname);	
      }
    }
    else {
      strcpy(seqname, tmp_seqname);
    }
    
    if ( nfields != 9 )
      error("The number of fields are not exactly 9\n%s", tr.str.s);

    // const char* source  = tr.str_field_at(1);
    const char* feature = tr.str_field_at(2);
    int32_t start       = tr.int_field_at(3);  // gtf beg is 1-inclusive
    int32_t end         = tr.int_field_at(4);  // gtf end is 1-inclusive
    const char* strand  = tr.str_field_at(6);
    const char* frame   = tr.str_field_at(7);
    const char* attr    = tr.str_field_at(8);

    if ( proteinCodingOnly ) {
      // search for gene_type "protein_coding"
      const char* pFind = strstr(attr, "gene_type \"protein_coding\"");
      if ( pFind == NULL ) continue;
    }    

    // tokenize the attrbute string
    int32_t  nattrs = 0;
    int32_t* attrFields;
    kstring_t kattr;
    kattr.l = kattr.m = strlen(attr);
    kattr.s = strdup(attr);
    attrFields = ksplit(&kattr, ';', &nattrs);
    --nattrs;
    for(int32_t i=1; i < nattrs; ++i) 
      ++attrFields[i];

    //bool fwdStrand = ( strand[0] == '-' ) ? false : true;

    // expected fields to exist
    // feature = gene        : gene_id, gene_name, gene_biotype
    // feature = transcript  : transcript_id, gene_id, trabscript_biotype
    // feature = exon        : transcript_id, exon_id, exon_number
    // feature = UTR         : transcript_id
    // feature = CDS         : transcript_id, exon_id, exon_number
    // feature = start_codon : transcript_id, exon_id, exon_number
    // feature = stop_codon  : transcript_id, exon_id, exon_number
    //notice("feature = %s, start = %d, end = %d", feature, start, end);

    if ( strcmp(feature, "five_prime_utr") == 0 ) {
      feature = "UTR";
    }
    else if ( strcmp(feature, "three_prime_utr") == 0 ) {
      feature = "UTR";      
    }
    
    std::string gid, tid, eid, name, type;
    if ( strcmp(feature,"gene") == 0 ) {
      for(int32_t i=0; i < nattrs; ++i) {
	if ( strncmp(&attr[attrFields[i]], "gene_id \"", 9) == 0 ) {
	  gid.assign(&attr[attrFields[i]+9], attrFields[i+1]-attrFields[i]-12);
	}
	else if ( strncmp(&attr[attrFields[i]], "gene_name \"", 11) == 0 ) {
	  name.assign(&attr[attrFields[i]+11], attrFields[i+1]-attrFields[i]-14);	  
	}
	else if ( strncmp(&attr[attrFields[i]], "gene_type \"", 11) == 0 ) {
	  type.assign(&attr[attrFields[i]+11], attrFields[i+1]-attrFields[i]-14);	  
	}
      }
      addGene(seqname, start, end, strand, gid, name, type);
    }
    else if ( strcmp(feature,"transcript") == 0 ) {
      for(int32_t i=0; i < nattrs; ++i) {
	if ( strncmp(&attr[attrFields[i]], "gene_id \"", 9) == 0 ) {
	  gid.assign(&attr[attrFields[i]+9], attrFields[i+1]-attrFields[i]-12);
	}
	else if ( strncmp(&attr[attrFields[i]], "transcript_id \"", 15) == 0 ) {
	  tid.assign(&attr[attrFields[i]+15], attrFields[i+1]-attrFields[i]-18);	  
	}
	else if ( strncmp(&attr[attrFields[i]], "transcript_type \"", 17) == 0 ) {
	  type.assign(&attr[attrFields[i]+17], attrFields[i+1]-attrFields[i]-20);	  
	}
      }
      //notice("%s %d %d %s %s %s %s\n%s",seqname, start, end, strand, gid.c_str(), tid.c_str(), type.c_str(), attr);
      addTranscript(seqname, start, end, strand, gid, tid, type);
    }
    else if ( strcmp(feature,"exon") == 0 ) {
      for(int32_t i=0; i < nattrs; ++i) {
	if ( strncmp(&attr[attrFields[i]], "exon_id \"", 9) == 0 ) {
	  eid.assign(&attr[attrFields[i]+9], attrFields[i+1]-attrFields[i]-12);
	}
	else if ( strncmp(&attr[attrFields[i]], "transcript_id \"", 15) == 0 ) {
	  tid.assign(&attr[attrFields[i]+15], attrFields[i+1]-attrFields[i]-18);	  
	}
      }
      addExon(seqname, start, end, strand, tid);
    }
    else if ( strcmp(feature,"UTR") == 0 ) {
      for(int32_t i=0; i < nattrs; ++i) {
	if ( strncmp(&attr[attrFields[i]], "transcript_id \"", 15) == 0 ) {
	  tid.assign(&attr[attrFields[i]+15], attrFields[i+1]-attrFields[i]-18);	  
	}
      }
      addUTR(seqname, start, end, strand, tid);
    }
    else if ( strcmp(feature,"CDS") == 0 ) {
      for(int32_t i=0; i < nattrs; ++i) {
	if ( strncmp(&attr[attrFields[i]], "transcript_id \"", 15) == 0 ) {
	  tid.assign(&attr[attrFields[i]+15], attrFields[i+1]-attrFields[i]-18);	  
	}
      }
      addCDS(seqname, start, end, strand, frame, tid);
    }
    else if ( strcmp(feature,"start_codon") == 0 ) {
      for(int32_t i=0; i < nattrs; ++i) {
	if ( strncmp(&attr[attrFields[i]], "transcript_id \"", 15) == 0 ) {
	  tid.assign(&attr[attrFields[i]+15], attrFields[i+1]-attrFields[i]-18);	  
	}
      }
      addStartCodon(seqname, start, end, strand, tid);
    }
    else if ( strcmp(feature,"stop_codon") == 0 ) {
      for(int32_t i=0; i < nattrs; ++i) {
	if ( strncmp(&attr[attrFields[i]], "transcript_id \"", 15) == 0 ) {
	  tid.assign(&attr[attrFields[i]+15], attrFields[i+1]-attrFields[i]-18);	  
	}
      }
      addStopCodon(seqname, start, end, strand, tid);
    }
    free(kattr.s);

    //notice("end of loop");    
  }
  notice("Finished reading %d lines from %s", line, filename);

  notice("Building interval trees for each contig..");
  for(gtf_chr_it_t it = mmap.begin(); it != mmap.end(); ++it) {
    notice("Processing contig %s", it->first.c_str());
    gtf_ivt_t::gtf_interval_vector itvs;
    for(gtf_elem_it_t jt = it->second.begin(); jt != it->second.end(); ++jt) {
      itvs.emplace_back(jt->first.beg1, jt->first.end0, jt->second);
    }
    //gtf_ivt_t newTree;
    //newTree = gtfIntervalTree<int32_t,gtfElement*>(itvs);
    //gtf_ivt_t::gtf_interval_vector& ritvs = itvs;
    chr2ivt[it->first] = gtf_ivt_t(std::move(itvs)); // std::move is necessary here
  }
}

gtf::~gtf() {
  // iterate across all elements, delete them
  notice("Deleting a GTF object");
  /*
  for(gtf_chr_it_t cit = mmap.begin(); cit != mmap.end(); ++cit) {
    gtf_chr_t& cmap = mit->second;
    
    for(gtf_elem_it_t eit = cmap.begin(); eit != cmap.end(); ++eit) {
      delete eit->second;
    }
  }
  */
  for(rewind(); !isend(); next()) {
    delete curElemIt->second;
  }
  
  notice("Finished deleting a GTF object");  
}

// create a gene
bool gtf::addGene(const char* seqname, int32_t start, int32_t end,
		  const char* strand, 
		  std::string& gid, std::string& gname, std::string& gtype) {
  // check if the gid already exists
  if ( gid2Gene.find(gid) != gid2Gene.end() ) {
    error("[E:%s:%d:%s] Duplicated gene ID %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, gid.c_str()); 
  }

  // figure out the strand info
  bool fwdStrand = true;
  if ( strand[0] == '+' )      { fwdStrand = true; }
  else if ( strand[0] == '-' ) { fwdStrand = false; }
  else error("[E:%s:%d:%s] Unknown strand info %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, strand);

  gtfGene* pGene = new gtfGene(seqname, start, end, fwdStrand, gid, gname, gtype);

  // add gene id into dictionary
  gid2Gene[gid] = pGene;
  // add the entry to the element map
  mmap[seqname].emplace(pGene->locus,pGene);
  //chr2ivt[seqname].push_back(gtf_iv_t(start, end, pGene));
  
  maxGeneLength = maxGeneLength > pGene->locus.length() ? maxGeneLength : pGene->locus.length();

  //notice("Added gene %s", gid.c_str());

  return true;
}

bool gtf::addTranscript(const char* seqname, int32_t start, int32_t end,
			const char* strand,
			std::string& gid, std::string& tid, std::string& ttype) {
  // check if the gid already exists
  auto git = gid2Gene.find(gid);
  if (git == gid2Gene.end() ) {
    error("[E:%s:%d:%s] Missing gene ID %s when inserting transcript %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, gid.c_str(), tid.c_str()); 
  }

  // check if the transcript ID was already observed
  if ( tid2Transcript.find(tid) != tid2Transcript.end() ) {
    error("[E:%s:%d:%s] Duplicated transcript ID %s", tid.c_str()); 
  }
  
  // create a new transcript
  gtfTranscript* pTranscript = new gtfTranscript(start, end, tid, ttype, git->second);
  auto ret = git->second->transcripts.emplace(pTranscript);
  if ( !ret.second )
    error("[E:%s:%d:%s] Duplicated element in transcript %s at %s:%d-%d", tid.c_str(), seqname, start, end);     

  // add transcript id into dictionary
  tid2Transcript[tid] = pTranscript;
  tid2Gene[tid] = git->second;
  
  // add the entry to the element map
  mmap[seqname].emplace( pTranscript->locus, pTranscript );
  //chr2ivt[seqname].push_back(gtf_iv_t(start, end, pTranscript));  
  maxTranscriptLength = maxTranscriptLength > pTranscript->locus.length() ? maxTranscriptLength : pTranscript->locus.length();

  //notice("Added transcript %s", tid.c_str());  

  return checkTranscriptSanity(tid, seqname, strand);        
}

bool gtf::addExon(const char* seqname, int32_t start, int32_t end,
		  const char* strand, 
		  std::string& tid) {
  // check if the gid already exists
  auto tit = tid2Transcript.find(tid);
  if (tit == tid2Transcript.end() ) {
    error("[E:%s:%d:%s] Missing transcript ID %s when inserting a new exon", __FILE__, __LINE__, __PRETTY_FUNCTION__, tid.c_str());
  }

  // add a new exon as a new gtfElement
  gtfElement* pExon = new gtfElement(start, end, "exon", tit->second);
  auto ret = tit->second->exons.emplace(pExon);
  if ( !ret.second )
    error("[E:%s:%d:%s] Duplicated element in transcript %s at %s:%d-%d", tid.c_str(), seqname, start, end);     

  // add the entry to the element map
  mmap[seqname].emplace( pExon->locus, pExon );
  //chr2ivt[seqname].push_back(gtf_iv_t(start, end, pExon));    
  maxExonLength = maxExonLength > pExon->locus.length() ? maxExonLength : pExon->locus.length();    

  return checkTranscriptSanity(tid, seqname, strand);      
}

bool gtf::addUTR(const char* seqname, int32_t start, int32_t end,
		  const char* strand, 
		  std::string& tid) {
  // check if the gid already exists
  auto tit = tid2Transcript.find(tid);
  if (tit == tid2Transcript.end() ) {
    error("[E:%s:%d:%s] Missing transcript ID %s when inserting a new UTR", __FILE__, __LINE__, __PRETTY_FUNCTION__, tid.c_str());
  }

  // add a new UTR as a new gtfElement
  gtfElement* pUTR = new gtfElement(start, end, "UTR", tit->second);
  auto ret = tit->second->UTRs.emplace(pUTR);
  if ( !ret.second )
    error("[E:%s:%d:%s] Duplicated element in transcript %s at %s:%d-%d", tid.c_str(), seqname, start, end);       

  // add the entry to the element map
  mmap[seqname].emplace( pUTR->locus, pUTR );
  //chr2ivt[seqname].push_back(gtf_iv_t(start, end, pExon));
  maxUTRLength = maxUTRLength > pUTR->locus.length() ? maxUTRLength : pUTR->locus.length();    

  return checkTranscriptSanity(tid, seqname, strand);      
}

bool gtf::addCDS(const char* seqname, int32_t start, int32_t end,
		 const char* strand, const char* frame,
		 std::string& tid)
{
  // check if the gid already exists
  auto tit = tid2Transcript.find(tid);
  if (tit == tid2Transcript.end() ) {
    error("[E:%s:%d:%s] Missing transcript ID %s when inserting a new CDS", __FILE__, __LINE__, __PRETTY_FUNCTION__, tid.c_str());
  }

  // add a new CDS as a new gtfElement
  gtfCDS* pCDS = new gtfCDS(start, end, frame, tit->second);
  auto ret = tit->second->CDSs.emplace(pCDS);
  if ( !ret.second )
    error("[E:%s:%d:%s] Duplicated element in transcript %s at %s:%d-%d", tid.c_str(), seqname, start, end);         

  // add the entry to the element map
  mmap[seqname].emplace( pCDS->locus, pCDS );
  //chr2ivt[seqname].push_back(gtf_iv_t(start, end, pCDS));  
  maxCDSLength = maxCDSLength > pCDS->locus.length() ? maxCDSLength : pCDS->locus.length();    

  return checkTranscriptSanity(tid, seqname, strand);      
}

bool gtf::addStartCodon(const char* seqname, int32_t start, int32_t end,
			const char* strand,
			std::string& tid)
{
  // check if the gid already exists
  auto tit = tid2Transcript.find(tid);
  if (tit == tid2Transcript.end() ) {
    error("[E:%s:%d:%s] Missing transcript ID %s when inserting a new CDS", __FILE__, __LINE__, __PRETTY_FUNCTION__, tid.c_str());
  }

  // add a new CDS as a new gtfElement
  gtfElement* pElement = new gtfElement(start, end, "start_codon", tit->second);  
  auto ret = tit->second->start_codons.emplace(pElement);
  if ( !ret.second )
    error("[E:%s:%d:%s] Duplicated element in transcript %s at %s:%d-%d", tid.c_str(), seqname, start, end);           

  // add the entry to the element map
  mmap[seqname].emplace( pElement->locus, pElement );
  //chr2ivt[seqname].push_back(gtf_iv_t(start, end, pElement));  
  maxStartCodonLength = maxStartCodonLength > pElement->locus.length() ? maxStartCodonLength : pElement->locus.length();    

  return checkTranscriptSanity(tid, seqname, strand);
}

bool gtf::addStopCodon(const char* seqname, int32_t start, int32_t end,
			const char* strand,
			std::string& tid)
{
  // check if the gid already exists
  auto tit = tid2Transcript.find(tid);
  if (tit == tid2Transcript.end() ) {
    error("[E:%s:%d:%s] Missing transcript ID %s when inserting a new CDS", __FILE__, __LINE__, __PRETTY_FUNCTION__, tid.c_str());
  }

  // add a new CDS as a new gtfElement
  gtfElement* pElement = new gtfElement(start, end, "stop_codon", tit->second);  
  auto ret = tit->second->stop_codons.emplace(pElement);
  if ( !ret.second )
    error("[E:%s:%d:%s] Duplicated element in transcript %s at %s:%d-%d", tid.c_str(), seqname, start, end);           

  // add the entry to the element map
  mmap[seqname].emplace( pElement->locus, pElement );
  //chr2ivt[seqname].push_back(gtf_iv_t(start, end, pElement));  
  maxStopCodonLength = maxStopCodonLength > pElement->locus.length() ? maxStopCodonLength : pElement->locus.length();    

  return checkTranscriptSanity(tid, seqname, strand);
}

bool gtf::checkTranscriptSanity(std::string& tid, const char* seqname, const char* strand) {
  // sanity check on seqname -- should match with seqname of gene
  auto git = tid2Gene.find(tid);
  if (git == tid2Gene.end() ) {
    error("[E:%s:%d:%s] Missing transcript ID %s when inserting a new CDS", __FILE__, __LINE__, __PRETTY_FUNCTION__, tid.c_str());
  }
  
  if ( git->second->seqname.compare(seqname) != 0 )
    error("[E:%s:%d:%s] seqname mismatches between gid %s (%s) and tid %s (%s)", __FILE__, __LINE__, __PRETTY_FUNCTION__, git->second->geneId.c_str(), git->second->seqname.c_str(), tid.c_str(), seqname);

  // sanity check on strand -- should match with strand of gene
  bool fwdStrand = true;
  if ( strand[0] == '+' )      { fwdStrand = true; }
  else if ( strand[0] == '-' ) { fwdStrand = false; }
  else error("[E:%s:%d:%s] Unknown strand info %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, strand);
  if ( git->second->fwdStrand != fwdStrand )
    error("[E:%s:%d:%s] strand mismatches between gid %s (%d) and tid %s (%d)", __FILE__, __LINE__, __PRETTY_FUNCTION__, git->second->geneId.c_str(), git->second->fwdStrand, tid.c_str(), fwdStrand);    

  return true;      
}

void gtf::rewind() {
  curChrIt = mmap.begin();
  curElemIt = curChrIt->second.begin();
}

bool gtf::isend() const {
  return curChrIt == mmap.end();
}

bool gtf::next() {
  ++curElemIt;
  if ( curElemIt == curChrIt->second.end() ) {
    while ( curChrIt != mmap.end() ) {
      ++curChrIt;
      curElemIt = curChrIt->second.begin();
      if ( curElemIt != curChrIt->second.end() )
	return true;
    }
    return false;
  }
  else
    return true;
}

int32_t gtf::findOverlappingElements(const char* seqname, int32_t start, int32_t end, std::set<gtfElement*>& results) {
  std::map<std::string, gtf_ivt_t>::iterator chr2ivt_it_t = chr2ivt.find(seqname);  
  if ( chr2ivt_it_t == chr2ivt.end() ) {
    notice("WARNING: no overlapping elements in %s", seqname);
    return 0;
  }
  gtf_ivt_t::gtf_interval_vector overlaps = chr2ivt_it_t->second.findOverlapping(start, end);
  for(int32_t i=0; i < (int32_t)overlaps.size(); ++i) {
    const gtf_ivt_t::gtf_interval& iv = overlaps[i];
    results.insert(iv.value);
    /*
    gtfElement* e = iv.value;
    gtfElement* root = e;
    while( root->parent != NULL )
      root = root->parent;
    gtfGene* rootGene = (gtfGene*)root;
    notice("Query = %s:%d-%d, Found %s:%d-%d, Type = %s, Gene ID = %s, Gene Name = %s, Gene Type = %s", seqname, start, end, seqname, iv.start, iv.stop, e->type.c_str(), rootGene->geneId.c_str(), rootGene->geneName.c_str(), rootGene->geneType.c_str());
    */
  }
  return 0;
}
