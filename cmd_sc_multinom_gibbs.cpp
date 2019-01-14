#include "cramore.h"
#include "tsv_reader.h"
#include "discrete_log_helper.h"

int32_t cmdScMultinomGibbs(int32_t argc, char** argv) {
  std::string inMatrix;
  std::string mtxf;
  std::string bcdf;
  std::string genef;    
  std::string outPrefix;
  //double doublet = 0;       // doublet proability
  int32_t burnin = 5;
  int32_t maxIter = 100;
  int32_t thin = 1;
  int32_t nClust = 0;       // Number of clusters required
  int32_t seed = 0;         // random seed
  int32_t geneThres = 1;
  double alpha = 1.0;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Required Options", NULL)
    LONG_STRING_PARAM("in",&inMatrix, "Input matrix int the format of R-compatible text matrix (can be gzipped)")
    LONG_STRING_PARAM("mtx",&mtxf, "Spare matrix representation (in .mtx format)")
    LONG_STRING_PARAM("bcd",&bcdf, "Barcode file used with --mtx option (in .tsv format)")
    LONG_STRING_PARAM("gene",&genef, "Barcode file used with --mtx option (in .tsv format)")            
    LONG_STRING_PARAM("out",&outPrefix, "Output file prefix")
    LONG_INT_PARAM("k",&nClust, "Number of clusters")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("burnin",&burnin,"Burn-in round")
    LONG_INT_PARAM("thin",  &thin,"Thinning interval")    
    LONG_INT_PARAM("gene-thres",&geneThres,"Burn-in round")    
    //LONG_DOUBLE_PARAM("doublet", &doublet, "Probability of being doublet")
    LONG_DOUBLE_PARAM("alpha",&alpha, "Pseudo-count per cell")
    //LONG_DOUBLE_PARAM("thres",&thresDiff, "Threshold of LLK difference to terminate the EM iteration")
    LONG_INT_PARAM("max-iter",&maxIter, "Number of maximum E-M iterations")
    LONG_INT_PARAM("seed",&seed, "Seed for random number generator (default uses clock)")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( nClust == 0 ) {
    error("[E:%s:%d %s] --k is a required parameter",__FILE__,__LINE__,__FUNCTION__);
  }

  if ( outPrefix.empty() ) {
    error("[E:%s:%d %s] --out is a required parameter",__FILE__,__LINE__,__FUNCTION__);
  }

  if ( inMatrix.empty() && ( mtxf.empty() || bcdf.empty() || genef.empty() ) ) {
    error("[E:%s:%d %s] --in or --mtx/--bcd/--gene is a required parameter",__FILE__,__LINE__,__FUNCTION__);
  }

  // open file to write
  htsFile* wf = hts_open((outPrefix+".pis").c_str(),"w");

  // read and parse the matrix elements
  int64_t nZero = 0;
  int64_t nSum = 0;
  int32_t nEmptyRows = 0;
  
  std::vector<std::string> hdrs;
  std::vector<std::string> genes;
  std::vector<int32_t*> R;
  std::vector<int64_t> rowSums;
  int64_t* colSums = NULL;  
  
  if ( wf == NULL )
    error("[E:%s:%d %s] Cannot open file %s for writing",__FILE__,__LINE__,__FUNCTION__, (outPrefix+".pis").c_str());

  // open input matrix
  if ( !inMatrix.empty() ) {
    htsFile* hp = hts_open(inMatrix.c_str(), "r");
    if ( hp == NULL )
      error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__,inMatrix.c_str());
    
    kstring_t str = {0,0,0};  
    int32_t lstr = 0;
    
    // read and parse header columns
    lstr = hts_getline(hp, KS_SEP_LINE, &str);
    if ( lstr < 0 )
      error("[E:%s:%d %s] Cannot find header line from %s",__FILE__,__LINE__,__FUNCTION__,inMatrix.c_str());
    
    int32_t nfields = 0;
    int32_t* fields = NULL;  
    fields = ksplit(&str, 0, &nfields);

    for(int32_t i=0; i < nfields; ++i) {
      hdrs.push_back(std::string(&str.s[fields[i]]));
    }
    
    notice("%d columns found in the header",(int32_t)hdrs.size());
    
    while( ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0 ) {
      fields = ksplit(&str, 0, &nfields);
      
      if ( R.empty() ) {
	notice("%d columns found in the second line",nfields);      
	if ( nfields == (int32_t)hdrs.size() ) { // remove one from the header
	  hdrs.erase(hdrs.begin());
	  notice("Ignoring the first column in the header");	
	}
	colSums = (int64_t*)calloc((int32_t)hdrs.size(),sizeof(int64_t));
      }
      else {
	if ( ( nfields != (int32_t)hdrs.size() + 1 ) && ( nfields != (int32_t)hdrs.size() + 0 ) )
	  error("[E:%s:%d %s] Inconsistent number of headers. Expected %d but observed %d",__FILE__,__LINE__,__FUNCTION__,(int32_t)hdrs.size()+1, nfields);
      }
      
      int32_t* cnts = (int32_t*)malloc(sizeof(int32_t)*(nfields-1));
      int64_t rowSum = 0;
      for(int32_t i=1; i < nfields; ++i) {
	cnts[i-1] = atoi(&str.s[fields[i]]);

	/*
	if ( fracSubsample < 1 ) {
	  int32_t tot = cnts[i-1];
	  int32_t sampled = 0;
	  for(int32_t j=0; j < tot; ++j) {
	    if ( (rand()+0.5) / (RAND_MAX+1.) < fracSubsample )
	      ++sampled;
	  }
	  cnts[i-1] = sampled;
	}
	*/
	
	if ( cnts[i-1] == 0 ) ++nZero;
	else {
	  rowSum += cnts[i-1];
	  colSums[i-1] += cnts[i-1];
	}
      }
      
      if ( rowSum < geneThres ) {
	if ( geneThres > 1 ) {
	  for(int32_t i=1; i < nfields; ++i) {
	    if ( cnts[i-1] == 0 ) --nZero;
	    colSums[i-1] -= cnts[i-1];
	  }
	}
	
	free(cnts);
	++nEmptyRows;
      }
      else {
	genes.push_back(std::string(&str.s[fields[0]]));      
	R.push_back(cnts);
	rowSums.push_back(rowSum);
	nSum += rowSum;
      }
    }
    hts_close(hp);
  }
  else {
    tsv_reader tr(bcdf.c_str());
    while(tr.read_line() > 0) {
      if ( tr.nfields > 0 )
	hdrs.push_back(tr.str_field_at(0));
    }

    tsv_reader tr2(genef.c_str());
    while(tr2.read_line() > 0) {
      if ( tr2.nfields > 0 )      
	genes.push_back(std::string(tr2.str_field_at(0)) + "_" + tr2.str_field_at(1) );
    }

    rowSums.resize(genes.size(),0);
    colSums = (int64_t*)calloc((int32_t)hdrs.size(),sizeof(int64_t));
    nSum = 0;
    nZero = (int64_t)genes.size() * (int64_t)hdrs.size();
    for(int32_t i=0; i < (int32_t)genes.size(); ++i) {
      R.push_back( (int32_t*) calloc(sizeof(int32_t), (int32_t)hdrs.size()) );
    }

    tsv_reader tr3(mtxf.c_str());
    tr3.read_line();
    tr3.read_line();
    tr3.read_line();
    if ( (int32_t)hdrs.size() != tr3.int_field_at(1) ) {
      error("Number of barcodes mismatch");
    }
    while(tr3.read_line() > 0) {
      if ( tr3.nfields > 0 )  {
	int32_t igene = tr3.int_field_at(0);
	int32_t ibcd =  tr3.int_field_at(1);
	int32_t cnt =  tr3.int_field_at(2);
	--nZero;
	//nSum += cnt;
	R[igene-1][ibcd-1] += cnt;
	rowSums[igene-1] += cnt;
	colSums[ibcd-1] += cnt;
      }
    }

    nEmptyRows = 0;
    
    for(int32_t i=0; i < (int32_t)genes.size(); ++i) {
      if ( rowSums[i] < geneThres ) {
	++nEmptyRows;
	if ( geneThres > 1 ) {
	  for(int32_t j=0; j < (int32_t)hdrs.size(); ++j) {
	    if ( R[i][j] > 0 ) {
	      ++nZero;
	      colSums[j] -= R[i][j];
	    }
	  }
	}
	free(R[i]);
	continue;
      }
      else if ( nEmptyRows > 0 ) {
	genes[i-nEmptyRows] = genes[i];
	R[i-nEmptyRows] = R[i];
	rowSums[i-nEmptyRows] = rowSums[i];
      }
      nSum += rowSums[i-nEmptyRows];
    }
    genes.resize(genes.size()-nEmptyRows);
    R.resize(genes.size());
    rowSums.resize(genes.size());
    
    notice("nEmptyRows = %d", nEmptyRows);
  }


  int32_t nRow = (int32_t)genes.size();
  int32_t nCol = (int32_t)hdrs.size(); // nCol is # of barcodes
  int64_t nCell = (int64_t)nRow * (int64_t)nCol;

  notice("Loaded a matrix with %d rows and %d columns after ignoring %d empty rows. Sparsity is %.5lg. Average of non-empty cells is %.5lg", nRow, nCol, nEmptyRows, (double)nZero/(double)(nCell+nEmptyRows*nCol), (double)nSum/(double)(nCell+nEmptyRows*nCol-nZero));

  if ( seed == 0 )
    srand(time(NULL));
  else
    srand(seed);  

  int32_t* zs =   (int32_t*) calloc ( nCol, sizeof(int32_t) );
  int32_t* knts = (int32_t*) calloc ( nClust, sizeof(int32_t) );
  int32_t* gsum = (int32_t*) calloc ( nClust * nRow, sizeof(int32_t) );
  int32_t* sumlg = (int32_t*) calloc ( nClust, sizeof(int32_t) );  
  
  double alpha_gi = alpha / nRow;
  double alpha_gk = alpha / nClust / nRow * nCol;
  double alpha_ik = alpha / nClust;  
  double alpha_i  = alpha;
  double alpha_k  = alpha / nClust * nCol;
  double alpha_a  = alpha * nCol;  
  discrete_log_helper dlh_gi(alpha_gi);  
  discrete_log_helper dlh_gk(alpha_gk);
  discrete_log_helper dlh_ik(alpha_ik);  
  discrete_log_helper dlh_i(alpha_i);
  discrete_log_helper dlh_k(alpha_k);
  discrete_log_helper dlh_a(alpha_a);  

  double* l_pi = (double*) calloc(nClust, sizeof(double));

  std::vector< std::vector<int32_t> > nzi(nCol);
  for(int32_t r=0; r < nRow; ++r) {
    for(int32_t c=0; c < nCol; ++c) {
      if ( R[r][c] > 0 ) nzi[c].push_back(r);
    }
  }

  // pick clusters randomly from the observations
  std::vector<int32_t> irand(nCol);
  for(int32_t c=0; c < nCol; ++c) irand[c] = c;
  for(int32_t k=0; k < nClust; ++k) {
    int32_t u = k + (int32_t)floor ( (rand() + 0.5)/(RAND_MAX + 1.0) * ( nCol - k ) );
    int32_t tmp = irand[u]; irand[u] = irand[k]; irand[k] = tmp;
  }

  for(int32_t c=0; c < nCol; ++c) {
    double max_pi = -1e300;
    for(int32_t k=0; k < nClust; ++k) {
      l_pi[k] = 0; // flat prior
      for(int32_t r=0; r < nRow; ++r) {
	l_pi[k] += ( (R[r][c] + alpha_gi) * dlh_gi.get(R[r][irand[k]]) );
      }
      l_pi[k] -= ( ( colSums[c] + alpha_i ) * dlh_i.get(colSums[irand[k]]) );
      
      if ( max_pi < l_pi[k] ) max_pi = l_pi[k];
    }
  
    double psum = 0;
    for(int32_t k=0; k < nClust; ++k) {
      psum += ( l_pi[k] = exp(l_pi[k] - max_pi) );
    }
    
    double u = (rand() + 0.5)/(RAND_MAX + 1.0) * psum;
    zs[c] = nClust-1;
    for(int32_t k=0; k < nClust; ++k) {
      if ( ( u -= l_pi[k] ) < 0 ) {
	zs[c] = k;
	break;
      }
    }
    if ( zs[c] >= nClust ) abort();
  }

  for(int32_t c=0; c < nCol-1; ++c) {
    ++knts[zs[c]];
  }

  notice("%d %d %d %d", knts[0], knts[1], knts[2], knts[3]);  

  for(int32_t c=0; c < nCol-1; ++c) {
    std::vector<int32_t>& vi = nzi[c];    
    for(int32_t r=0; r < (int32_t)vi.size(); ++r) {
      gsum[zs[c] * nRow + vi[r]] += R[vi[r]][c];
    }
  }

  for(int32_t k=0; k < nClust; ++k) {
    sumlg[k] = 0;
    for(int32_t r=0; r < nRow; ++r) {
      sumlg[k] += dlh_gk.get(gsum[k * nRow + r]);
    }
  }

  notice("Finished initializing the Gibbs sampler");  

  for(int32_t iter=0; iter < maxIter; ++iter) {
    for(int32_t k=0; k < nClust; ++k) {
      sumlg[k] = 0;
      for(int32_t r=0; r < nRow; ++r) {
	sumlg[k] += dlh_gk.get(gsum[k * nRow + r]);
      }
    }
    
    for(int32_t c=0; c < nCol; ++c) {
      int32_t o = (c == 0 ? nCol-1 : c-1);
      ++knts[zs[o]];
      --knts[zs[c]];

      std::vector<int32_t>& vi = nzi[o];
      for(int32_t r=0; r < (int32_t)vi.size(); ++r) {
	gsum[zs[o]*nRow+vi[r]] += R[vi[r]][o];
	sumlg[zs[o]] += ( dlh_gk.get(gsum[zs[o]*nRow+vi[r]]) - dlh_gk.get( gsum[zs[o]*nRow+vi[r]] - R[vi[r]][o] ) );
      }

      vi = nzi[c];
      for(int32_t r=0; r < (int32_t)vi.size(); ++r) {
	if ( ( gsum[zs[c]*nRow+vi[r]] -= R[vi[r]][c] ) < 0 )
	  error("r = %d, vi[r] = %d, c = %d, zs[c] = %d, gsum = (%d %d %d %d), R = %d, rowSums = %d, zs[o] = %d, knts = (%d %d %d %d)", r, vi[r], c, zs[c], gsum[vi[r]], gsum[nRow+vi[r]], gsum[2*nRow+vi[r]], gsum[3*nRow+vi[r]], R[vi[r]][c], rowSums[vi[r]], zs[o], knts[0], knts[1], knts[2], knts[3]); 
	sumlg[zs[c]] += ( dlh_gk.get(gsum[zs[c]*nRow+vi[r]]) - dlh_gk.get( gsum[zs[c]*nRow+vi[r]] + R[vi[r]][c] ) );
      }

      double max_pi = -1e300;
      for(int32_t k=0; k < nClust; ++k) {
	l_pi[k] = dlh_k.get(knts[k]); // prior based on current state
	for(int32_t r=0; r < (int32_t)vi.size(); ++r) {
	  l_pi[k] += ( R[vi[r]][c] * dlh_gk.get(gsum[k*nRow+vi[r]]) );
	}
	l_pi[k] += ( alpha_i * sumlg[k] - ( alpha_i + colSums[c] ) * dlh_k.get(knts[k]) );
	if ( max_pi < l_pi[k] ) max_pi = l_pi[k];
      }

      double psum = 0;
      for(int32_t k=0; k < nClust; ++k) {
	psum += ( l_pi[k] = exp(l_pi[k] - max_pi) );
      }

      int32_t mknt = knts[0];
      int32_t imknt = 0;
      for(int32_t k=1; k < nClust; ++k) {
	if ( knts[k] < mknt ) {
	  mknt = knts[k];
	  imknt = k;
	}
      }

      if ( mknt < 10 ) {
	zs[c] = imknt;
      }
      else {
	double u = (rand() + 0.5)/(RAND_MAX + 1.0) * psum;
	zs[c] = nClust-1;
	for(int32_t k=0; k < nClust; ++k) {
	  if ( ( u -= l_pi[k] ) < 0 ) {
	    zs[c] = k;
	    k = nClust; //break;
	  }
	}
	if ( zs[c] >= nClust ) abort();
      }

      //if ( knts[zs[c]] < 10 )
      //notice("bar %d %d %d", c, zs[c], knts[zs[c]]);      

      if ( c % 1000 == 0 ) notice("%d %d %d %d %d", c, knts[0], knts[1], knts[2], knts[3]);
    }

    if ( ( iter > burnin ) && ( iter % thin == 0 ) ) {
      notice("\nIter = %d", iter);
      for(int32_t k=0; k < nClust; ++k) {
	notice("pis[%d] = %.3lg", k, exp(dlh_ik.get(knts[k]) - dlh_i.get(nCol-1)));
      }
    }
    else if ( iter % thin == 0 ) {
      notice("Finished burn-in iteration %d", iter);      
    }
  }

  free(l_pi);
  free(zs);
  free(knts);
  free(gsum);  
  
  return 0;
}

