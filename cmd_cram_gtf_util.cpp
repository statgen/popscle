#include "cramore.h"
#include "gtf.h"

int32_t cmdGtfUtil(int32_t argc, char** argv) {
  std::string inGTF;
  std::string outPrefix;
  bool summaryFlag;
  bool genePredFlag;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input GTF", NULL)
    LONG_STRING_PARAM("gtf",&inGTF, "Input GTF file. Plain, gzipped, or bgzipped formats are  allowed")

    LONG_PARAM_GROUP("Options for output file", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Output file name")
    
    LONG_PARAM_GROUP("Analysis Options", NULL)
    LONG_PARAM("summary",&summaryFlag, "Provide a concise summary of the GTF file")
    LONG_PARAM("gene-pred",&genePredFlag, "Convert the GTF file in UCSC genePreds format")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( inGTF.empty() || outPrefix.empty() )
    error("[E:%s:%d %s] --gtf or --out parameter is missing",__FILE__,__LINE__,__PRETTY_FUNCTION__);

  gtf myGTF(inGTF.c_str());

  if ( summaryFlag ) {
    notice("Max gene        length = %d", myGTF.maxGeneLength);
    notice("Max transcript  length = %d", myGTF.maxTranscriptLength);
    notice("Max exon        length = %d", myGTF.maxExonLength);
    notice("Max UTR         length = %d", myGTF.maxUTRLength);
    notice("Max CDS         length = %d", myGTF.maxCDSLength);
    notice("Max start_codon length = %d", myGTF.maxStartCodonLength);      
    notice("Max stop_codon  length = %d", myGTF.maxStopCodonLength);
  }

  if ( genePredFlag ) {
    htsFile* wf;
    if ( outPrefix.compare( outPrefix.size()-3, 3, ".gz" ) == 0 )
      wf = hts_open(outPrefix.c_str(), "wg");
    else
      wf = hts_open(outPrefix.c_str(), "w");      
    
    notice("Writing GTF file in genePred format to %s", outPrefix.c_str());
    int32_t ntr = 0;
    for( myGTF.rewind(); !myGTF.isend(); myGTF.next() ) {
      gtfElement* e = myGTF.curElemIt->second;
      if ( e->type == "transcript" ) {  // focus only on transcripts
	gtfTranscript* t = (gtfTranscript*) e;
	gtfGene* g = (gtfGene*)t->parent;
	hprintf(wf, "%s", t->transcriptId.c_str());      
	hprintf(wf, "\t%s", g->seqname.c_str());
	if ( g->fwdStrand ) hprintf(wf, "\t+");
	else                hprintf(wf, "\t-");
	hprintf(wf,"\t%d\t%d", t->locus.beg1-1, t->locus.end0);
	
	// check if CDS exists
	if ( t->CDSs.empty() ) {
	  hprintf(wf,"\t%d\t%d", t->locus.end0, t->locus.end0);	
	}
	else {
	  hprintf(wf, "\t%d", (*t->CDSs.begin())->locus.beg1-1);
	  hprintf(wf, "\t%d", (*t->CDSs.rbegin())->locus.end0);	
	}
	
	if ( t->exons.empty() )
	  error("FATAL ERROR: Empty exon count");
	hprintf(wf, "\t%u", t->exons.size());
	
	// print out exonStarts
	std::set<gtfElement*, gtfComp>::iterator eit = t->exons.begin();
	hprintf(wf, "\t");
	while( eit != t->exons.end() ) {
	  hprintf(wf, "%d,", (*eit)->locus.beg1-1);
	  ++eit;
	}
	
	// print out exonEnds
	eit = t->exons.begin();
	hprintf(wf, "\t");      
	while( eit != t->exons.end() ) {
	  hprintf(wf, "%d,", (*eit)->locus.end0);
	  ++eit;
	}
	
	hprintf(wf, "\t0"); // print score, which is always zero
	hprintf(wf, "\t%s", g->geneName.c_str());
	
	if ( t->CDSs.empty() ) {
	  hprintf(wf, "\tnone\tnone\t");
	}
	else {
	  if ( g->fwdStrand ) { 
	    if ( t->start_codons.empty() ) hprintf(wf, "\tincmpl");
	    else hprintf(wf, "\tcmpl");
	    if ( t->stop_codons.empty() ) hprintf(wf, "\tincmpl");
	    else hprintf(wf, "\tcmpl");	    
	  }
	  else {
	    if ( t->stop_codons.empty() ) hprintf(wf, "\tincmpl");
	    else hprintf(wf, "\tcmpl");
	    if ( t->start_codons.empty() ) hprintf(wf, "\tincmpl");
	    else hprintf(wf, "\tcmpl");	    	    
	  }
	}
	
	if ( t->CDSs.empty() ) {
	  eit = t->exons.begin();
	  while( eit != t->exons.end() ) {
	    hprintf(wf, "-1,");
	    ++eit;
	  }	
	}
	else {
	  std::set<gtfElement*,gtfComp>::iterator exonIt = t->exons.begin();
	  std::set<gtfCDS*,gtfComp>::iterator cdsIt  = t->CDSs.begin();
	  while( exonIt != t->exons.end() ) {
	    if ( cdsIt == t->CDSs.end() ) hprintf(wf, "-1,");
	    else if ( (*cdsIt)->locus.overlaps((*exonIt)->locus) ) {
	      hprintf(wf, "%d,", (*cdsIt)->frame);
	      ++cdsIt;
	    }
	    else hprintf(wf, "-1,");
	    ++exonIt;
	  }
	}
	
	hprintf(wf, "\t%s\n", g->geneId.c_str());
	
	++ntr;
	
	if ( ntr % 10000 == 0 )
	  notice("Writing %d transcripts...", ntr);
      }
    }
    hts_close(wf);
    notice("Finished writing %d transcripts in genePred format", ntr);
  }
  
  return 0;
}
