#include "cramore.h"
// map a 12-letter STAMPs accounting for errors and trimming
// when length is 12, it is using 24-bits in 2bit space
// here are the information collected
// 1. 24bit maps counting the number of reads (randomly resolve Ns)
// 2. Sort the maps based on the counts and generate raw count distribution

int32_t cmdScMapSTAMPs(int32_t argc, char** argv) {
  std::string inFastQ;
  std::string outPrefix;  
  int32_t bcLen = 12;
  int32_t bcMinLen = 9;
  int32_t umiLen = 8;
  int32_t trailLen = 3;
  
  paramList pl;
  
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input FASTQ Options", NULL)
    LONG_STRING_PARAM("fq",&inFastQ, "Input FASTQ file (plain-text or gzipped)")
    
    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("bc-len",&bcLen, "Barcode length")
    LONG_INT_PARAM("bc-min-len",&bcMinLen, "Minimum Barcode length retained during the trimming")
    LONG_INT_PARAM("umi-len",&umiLen, "UMI length")
    LONG_INT_PARAM("trail-len",&trailLen, "Length of trailing bases")
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outPrefix, "Output prefix")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  htsFile* hp = hts_open(inFastQ.c_str(), "r");
  if ( hp == NULL )
    error("Cannot open file %s for reading",inFastQ.c_str());
  
  // read FASTQ files
  notice("Scanning the FASTQ file first time to construct the barcode map");
  kstring_t lines[4] = { {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0} };
  int32_t lens[4];

  std::vector<uint32_t> nBarcodes(bcLen+1,0);
  std::vector<uint32_t*> cntBarcodes(bcLen+1,NULL);
  std::vector<uint32_t> maxCnts(bcLen+1,0);  

  for(int32_t i=bcMinLen; i <= bcLen; ++i) {
    nBarcodes[i] =  (1 << (i*2));
    cntBarcodes[i] = (uint32_t*) calloc(nBarcodes[i], sizeof(uint32_t) );
  }

  uint32_t n2i[256] = {0};
  //memset(n2i, 0, sizeof(uint32_t)*256);
  n2i[(int32_t)'A'] = 0; n2i[(int32_t)'a'] = 0;
  n2i[(int32_t)'C'] = 1; n2i[(int32_t)'c'] = 1;
  n2i[(int32_t)'G'] = 2; n2i[(int32_t)'g'] = 2;
  n2i[(int32_t)'T'] = 3; n2i[(int32_t)'t'] = 3;    

  int32_t nRead = 0;
  
  while( ( lens[0] = hts_getline(hp, KS_SEP_LINE, &lines[0]) ) > 0 ) {
    if ( nRead % 10000000 == 0 )
      notice("Reading %d FASTQ Records", nRead);
           
    for(int32_t i=1; i < 4; ++i) {
      if ( ( lens[i] = hts_getline(hp, KS_SEP_LINE, &lines[i]) ) <= 0 ) {
	error("FASTQ number of lines are %d mod 4", i);
      }
    }

    // sequence is stored in lines[1]
    uint32_t bc = 0;
    char* s = lines[1].s;
    if ( lines[1].l < (uint32_t)bcLen )
      error("lines[1].l = %d < bcLen = %d",lines[1].l,bcLen);
    
    for(int32_t j=0; j < bcLen; ++j) {
      bc = ( (bc << 2) + ( (s[j] == 'N') ? (rand() % 4) : n2i[(int32_t)s[j]] ) );
      if ( j+1 >= bcMinLen ) {
	++(cntBarcodes[j+1][bc]);
	if ( maxCnts[j+1] < cntBarcodes[j+1][bc] )
	  maxCnts[j+1] = cntBarcodes[j+1][bc];
      }
    }

    ++nRead;
  }

  notice("Loaded total of %d reads", nRead);

  // calculate the summary statistics of barcode counts
  std::vector< std::vector<uint32_t> > bcHist(bcLen+1);
  for(int32_t i=bcMinLen; i <= bcLen; ++i) {
    bcHist[i].resize(maxCnts[i]+1);
    for(uint32_t j = 0; j < nBarcodes[i]; ++j)
      ++(bcHist[i][cntBarcodes[i][j]]);

    notice("Writing raw barcode histogram for length %d - maxCnt is %d",i, maxCnts[i]);

    char buf[255];
    sprintf(buf,"%d",i);

    htsFile* wf = hts_open((outPrefix+".raw."+buf+".hist").c_str(),"w");
    uint32_t sum = 0;
    for(int32_t j=maxCnts[i]; j >= 0; --j) {
      if ( bcHist[i][j] > 0 ) {
	sum += bcHist[i][j];
	hprintf(wf, "%d\t%u\t%u\n", j, bcHist[i][j], sum);
	//fprintf(stderr, "%u\t%u\n", j, bcHist[i][j]);	
      }
    }
    hts_close(wf);
  }

  for(int32_t i=bcMinLen; i <= bcLen; ++i) {
    free(cntBarcodes[i]);
  }

  return 0;
}
