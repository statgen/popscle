#include "cramore.h"
#include "tsv_reader.h"

int32_t cmdScMultinomEM(int32_t argc, char** argv) {
  std::string inMatrix;
  std::string mtxf;
  std::string bcdf;
  std::string genef;    
  std::string outPrefix;
  double doublet = 0;       // doublet proability
  double alpha = 0.5;       // pseudo-count per cell / gene
  double thresDiff = 1e-10; // threshold to stop EM iteration
  int32_t maxIter = 100;    // maximum number of EM iteration
  int32_t nClust = 0;       // Number of clusters required
  int32_t seed = 0;         // random seed
  int32_t nCollapseGenes = 0; // collapse genes into a specific number
  double fracSubsample = 1;   // fraction of samples to thin the data
  int32_t geneThres = 1;

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
    LONG_INT_PARAM("gene-thres",&geneThres,"Threshold for per gene count")
    LONG_DOUBLE_PARAM("doublet", &doublet, "Probability of being doublet")
    LONG_DOUBLE_PARAM("alpha",&alpha, "Pseudo-count per cell")
    LONG_DOUBLE_PARAM("thres",&thresDiff, "Threshold of LLK difference to terminate the EM iteration")
    LONG_INT_PARAM("max-iter",&maxIter, "Number of maximum E-M iterations")
    LONG_INT_PARAM("collapse-genes",&nCollapseGenes,"Number of genes to be collapsed into to reduce parameter space")
    LONG_INT_PARAM("seed",&seed, "Seed for random number generator (default uses clock)")
    LONG_DOUBLE_PARAM("frac-subsample",&fracSubsample, "Fraction of samples to thin the data")
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
	
	if ( fracSubsample < 1 ) {
	  int32_t tot = cnts[i-1];
	  int32_t sampled = 0;
	  for(int32_t j=0; j < tot; ++j) {
	    if ( (rand()+0.5) / (RAND_MAX+1.) < fracSubsample )
	      ++sampled;
	  }
	  cnts[i-1] = sampled;
	}
	
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
    /*
    for(int32_t i=0; i < (int32_t)genes.size(); ++i) {
      if ( rowSums[i] == 0 )
	++nEmptyRows;
	}*/

    //std::vector<std::string> genes;
    //std::vector<int32_t*> R;
    //std::vector<int64_t> rowSums;
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

  if ( nCollapseGenes > 0 ) {
    std::vector< std::vector<int32_t> > group2Gene( nCollapseGenes );
    std::vector< int32_t > gene2Group( nRow, 0 );

    for(int32_t i=0; i < nRow; ++i) {
      int32_t g = (rand() % nCollapseGenes);
      group2Gene[g].push_back(i);
    }

    nEmptyRows = 0;
		
    for(int32_t i=nCollapseGenes-1; i >= 0; --i) {
      if ( group2Gene[i].empty() ) {
	++nEmptyRows;
	group2Gene.erase(group2Gene.begin() + i);
      }
      else {
	for(int32_t j=0; j < (int32_t)group2Gene[i].size(); ++j) {
	  gene2Group[group2Gene[i][j]] = i;
	}
      }
    }

    std::vector<std::string> newGenes(nCollapseGenes-nEmptyRows);
    std::vector<int64_t> newRowSums(nCollapseGenes-nEmptyRows, 0);
    std::vector<int32_t*> newR(nCollapseGenes-nEmptyRows, NULL);

    for(int32_t i=0; i < nRow; ++i) {
      int32_t g = gene2Group[i];
      if ( newGenes[g].empty() ) {
	newGenes[g] = genes[i];
	newR[g] = (int32_t*) calloc(sizeof(int32_t), nCol);
      }
      else {
	newGenes[g] += ",";
	newGenes[g] += genes[i];
      }
      newRowSums[g] += rowSums[i];
      for(int32_t j=0; j < nCol; ++j) {
	newR[g][j] += R[i][j];
      }
      free(R[i]);
    }

    genes = newGenes;
    rowSums = newRowSums;
    R = newR;
    nRow = nCollapseGenes-nEmptyRows;
    nCell = (int64_t)nRow * (int64_t)nCol;    

    notice("Collapsed the matrix with %d rows and %d columns after ignoring %d additional empty rows created during the random collpaing procedure", nRow, nCol, nEmptyRows);
  }

  // calculate the global proportion matrix
  double* p0 = new double[nRow];
  for(int32_t i=0; i < nRow; ++i) {
    p0[i] = (double)rowSums[i]/(double)nSum;
  }

  // create multiple copies of parameters for simultaneous EM
  int32_t nPair = ( doublet > 0 ? nClust * (nClust-1) / 2 : 0 );
  int32_t hasDoublet = ( doublet > 0 ? 1 : 0);
  double log_doublet = doublet > 0 ? log(doublet) : 0;
  double log_singlet = doublet > 0 ? log(1-doublet) : 0;  
  //double* pis = (double*)calloc( (nClust + hasDoublet),sizeof(double));
  double* pis = (double*)calloc( nClust,sizeof(double));    
  double* Ps = (double*)calloc( nClust * nRow, sizeof(double));
  double* Ls = (double*)calloc( nClust * nRow, sizeof(double));  
  double* Zs = (double*)calloc( (nClust + nPair) * nCol,sizeof(double));
  double llk = 0, llk0 = 0;

  if ( seed == 0 )
    srand(time(NULL));
  else
    srand(seed);

  // randomize class assignments
  for(int32_t c=0; c < nCol; ++c) {            // for each barcode, 
    double* z = &Zs[c * ( nClust + nPair ) ];    
    double u = (rand() + 0.5) / (RAND_MAX + 1.0);
    double u2 = (rand() + 0.5) / (RAND_MAX + 1.0);    
    if ( u < doublet ) {
      //z[(int32_t)(floor((rand()+0.5)/(RAND_MAX+1.)*nPair)) + nClust] = 1.;
      z[ (int32_t)( u2 * (double)nPair) + nClust ] = 1.0;
      
    }
    else {
      //z[(int32_t)(floor((rand()+0.5)/(RAND_MAX+1.)*nClust))] = 1.;
      z[ (int32_t)(u2 * (double)nClust) ] = 1.0;      
    }
  }

  //notice("foo");

  // run EM iteration
  for(int32_t iter=0; iter < maxIter; ++iter) {
    // At this stage..
    // Z   : nCol x (nClust+nPair) - posterior probability of each barcode being assigned to each cluster (or pair)
    // pis : prior probability of each subclasses.    
    // Ps  : nRow x nClust        -- probability of each gene-cluster pair
       
    // M-step for pi
    memset(pis, 0, nClust * sizeof(double));
    
    for(int32_t c=0; c < nCol; ++c) {
      double* z = &Zs[c * (nClust+nPair)];    
      for(int32_t k=0; k < nClust; ++k) {
	pis[k] += z[k];  // pi_k = \sum_c Pr(z_c = k)
      }
      if ( hasDoublet > 0 ) {
	//for(int32_t k=0; k < nPair; ++k) {      
	//pis[nClust] += z[nClust+k];
	//}
	for(int32_t k1=1; k1 < nClust; ++k1) {
	  for(int32_t k2=0; k2 < k1; ++k2) {
	    pis[k1] += ( z[k1*(k1-1)/2+k2 + nClust] / 2 );
	    pis[k2] += ( z[k1*(k1-1)/2+k2 + nClust] / 2 );	    
	  }
	}
      }
    }


    //for(int32_t k=0; k < nClust+hasDoublet; ++k) {
    //  notice("k=%d\tu_pi=%lg",k,pis[k]);
    //}    

    double sum = 0;
    double pi2sum = 0;
    //for(int32_t k=0; k < nClust+hasDoublet; ++k) {
    for(int32_t k=0; k < nClust; ++k) {      
      sum += pis[k];
    }
    //for(int32_t k=0; k < nClust+hasDoublet; ++k) {
    for(int32_t k=0; k < nClust; ++k) {      
      pis[k] = log( pis[k]/sum );
    }

    if ( hasDoublet > 0 ) {
      for(int32_t k1=1; k1 < nClust; ++k1) {
	for(int32_t k2=0; k2 < k1; ++k2) {
	  pi2sum += exp(pis[k1] + pis[k2]);
	}
      }      
    }

    //notice("pi2sum = %lg, pis = (%lg, %lg, %lg, %lg, %lg, %lg)", pi2sum, exp(pis[0]), exp(pis[1]), exp(pis[2]), exp(pis[3]), exp(pis[4]), exp(pis[5]));
    pi2sum = log(pi2sum);

    //for(int32_t k=0; k < nCxR; ++k) {
    //  notice("k=%d\tpi=%lg",k,exp(pis[k]));
    //}

    //notice("iter = %d", iter);    
    
    // M-step for P (without normalization)
    for(int32_t g=0; g < nRow; ++g) {
      double* p = &Ps[g * nClust];
      
      for(int32_t k=0; k < nClust; ++k)
	p[k] = 0;
      
      for(int32_t c=0; c < nCol; ++c) {
	double* z = &Zs[c * (nClust+nPair)];
	double r = R[g][c] + alpha;
	//double r = R[g][c] + p0[g]*alpha;	
	//double r = R[g][c] + alpha/nRow;
	for(int32_t k=0; k < nClust; ++k) {
	  //double t = z[k] * r;
	  p[k] += (z[k] * r); // not normalized   \Pr(x_g|z_c=k) \propt \sum_c R_gc Pr(z_c=k)
	}
	if ( hasDoublet > 0 ) {
	  for(int32_t k1=1; k1 < nClust; ++k1) {
	    for(int32_t k2=0; k2 < k1; ++k2) {
	      p[k1] += (z[nClust + k1*(k1-1)/2 + k2] * r / 2);
	      p[k2] += (z[nClust + k1*(k1-1)/2 + k2] * r / 2);	      
	    }
	  }
	}
      }
    }

    //notice("goo");

    // normalize P
    for(int32_t k=0; k < nClust; ++k) {
      double sumP = 0;
      for(int32_t g=0; g < nRow; ++g) {
	sumP += Ps[g*nClust + k];
      }
      
      for(int32_t g=0; g < nRow; ++g) {
	Ps[g*nClust + k] /= sumP;
      }
    }
    

    // transform p into logp
    for(int32_t g=0; g < nRow; ++g) {
      double* p = &(Ps[g * nClust]);
      double* l = &(Ls[g * nClust]);      
      for(int32_t k=0; k < nClust; ++k) {
	l[k] = log(p[k]); // + pis[k];  // pi*P in log-scale
      }
    }


    for(int32_t c=0; c < nCol; ++c) {
      double* z = &(Zs[c*(nClust+nPair)]);      
      for(int32_t k=0; k < nClust; ++k) {
	z[k] = log_singlet + pis[k]; // probability to belong k-th cluster
      }
      if ( hasDoublet > 0 ) {
	for(int32_t k1=1; k1 < nClust; ++k1) {
	  for(int32_t k2=0; k2 < k1; ++k2) {
	    // (1-doublet) * pi(k1) * pi(k2) * 2 / Normalize
	    // z[nClust + k1*(k1-1)/2+k2] = pis[nClust] + pis[k1] + pis[k2] - pi2sum;
	    z[nClust + k1*(k1-1)/2+k2] = log_doublet + pis[k1] + pis[k2] - pi2sum;	    
	  }
	}	
      }
    }

    //for(int32_t k=0; k < nClust+hasDoublet; ++k)
    for(int32_t k=0; k < nClust; ++k)      
      notice("pis[%d] = %lg", k, exp(pis[k]));

    if ( hasDoublet > 0 ) {
      for(int32_t k1=1; k1 < nClust; ++k1) {
	for(int32_t k2=0; k2 < k1; ++k2) {
	  //notice("pis[%d,%d] = %lg", k2, k1, exp(pis[nClust] + pis[k1] + pis[k2] - pi2sum));
	  notice("pis[%d,%d] = %lg", k2, k1, exp(log_doublet + pis[k1] + pis[k2] - pi2sum));	  
	}
      }
    }
    
    // E-step for Z : t(R) %*% logP
    for(int32_t g=0; g < nRow; ++g) {
      double* p = &(Ps[g * nClust]);
      double* l = &(Ls[g * nClust]);      
      for(int32_t c=0; c < nCol; ++c) {
	double r = R[g][c] + alpha;
	//double r = R[g][c] + p0[g]*alpha; 	
	//int32_t r = R[g][c];
	double* z = &(Zs[c*(nClust+nPair)]);
	for(int32_t k=0; k < nClust; ++k) {
	  z[k] += (r*l[k]);  // \log Pr(z_c = k | x) \propt \sum_g [ \log \Pr(R_gc|z_c=k) ] + \pi_k
	}
	if ( hasDoublet > 0 ) {
	  for(int32_t k1=1; k1 < nClust; ++k1) {
	    for(int32_t k2=0; k2 < k1; ++k2) {
	      //double avgp = (p[k1] + p[k2]) / 2; //( ( p[k1] > p[k2] ) ? ( p[k1] + log(1.0 + exp(p[k2]-p[k1])) ) : ( p[k2] + log(1.0 + exp(p[k1]-p[k2])) ) ) - log(2.0);
	      z[nClust + k1*(k1-1)/2 + k2] += (r * log( (p[k1]+p[k2]) / 2.0));
	    }
	  }	
	}
      }
    }

    if ( iter == 0 ) llk0 = -1e300;
    else llk0 = llk;
    llk = 0;

    for(int32_t c=0; c < nCol; ++c) {
      double* z = &(Zs[c*(nClust+nPair)]);      
      double maxZ = z[0];
      for(int32_t k=1; k < nClust+nPair; ++k) {
	if ( maxZ < z[k])
	  maxZ = z[k];
      }

      double sumZ = 0;
      for(int32_t k=0; k < nClust+nPair; ++k) {
	double zdiff = z[k] - maxZ;	  
	sumZ += (z[k] = exp(zdiff));
      }

      llk += (maxZ + log(sumZ));
      if ( std::isnan(llk) ) {
	notice("maxZ = %lg, sumZ = %lg, llk = %lg, z[0] = %lg", maxZ, sumZ, llk, z[0]);
	abort();
      }
    }

    double maxDiff = -1e300;
    notice("Iter:%d\tLLK=%.5lf\tDiff=%.5lg", iter, llk, llk-llk0); //, exp(pis[0]), exp(pis[nRestarts]));
    if ( maxDiff < llk-llk0 ) {
      maxDiff = llk-llk0;
    }

    if ( maxDiff < thresDiff ) {
      notice("All LLK differences are less than %.5lg < %.5lg", maxDiff, thresDiff);
      break;
    }

    if ( iter + 1 == maxIter )
      notice("Reached maximum iteration %d",maxIter);


    /*
    for(int32_t c=0; c < 10; ++c) {
      double* z = &Zs[c * (nClust+nPair)];    
      for(int32_t k=0; k < nClust+nPair; ++k) 
	printf("%.3lg ", z[k]);
      printf("\n");
    }
    */
  }

  // transform P to linear scale
  /*
  for(int32_t k=0; k < nClust; ++k) {
    double maxP = Ps[k];
    for(int32_t g=1; g < nRow; ++g) {
      if ( maxP < Ps[g*nClust + k] )
	maxP = Ps[g*nClust + k];
    }

    double sumP = 0;
    for(int32_t g=0; g < nRow; ++g) {
      sumP += (Ps[g*nClust + k] = exp(Ps[g*nClust + k] - maxP));
    }

    for(int32_t g=0; g < nRow; ++g) {
      Ps[g*nClust + k] /= sumP;
    }
  }
  */

  for(int32_t k=0; k < nClust + hasDoublet; ++k) 
    hprintf(wf, "%g\n",exp(pis[k]));
  hts_close(wf);

  wf = hts_open((outPrefix+".Ps").c_str(),"w");
  for(int32_t g=0; g < nRow; ++g) {
    hprintf(wf, "%s",genes[g].c_str());    
    for(int32_t k=0; k < nClust; ++k) {
      hprintf(wf, "\t%.5lg",Ps[g*nClust + k]); 
    }
    hprintf(wf, "\n");
  }
  hts_close(wf);

  wf = hts_open((outPrefix+".Zs").c_str(),"w");
  for(int32_t c=0; c < nCol; ++c) {
    double sumZ = 0;
    for(int32_t k=0; k < nClust+nPair; ++k)
      sumZ += Zs[c*(nClust+nPair) + k];

    int32_t iBest = 0;
    for(int32_t k=1; k < nClust+nPair; ++k) {
      if ( Zs[c*(nClust+nPair) + iBest] < Zs[c*(nClust+nPair) + k] )
	iBest = k;
    }

    if ( iBest < nClust )
      hprintf(wf, "%s\t%d\t%d",hdrs[c].c_str(), colSums[c], iBest+1);
    else {
      int32_t k1, k2;
      for(k1=1; k1 < nClust; ++k1) {
	for(k2=0; k2 < k1; ++k2) {
	  if ( k1*(k1-1)/2 + k2 == iBest - nClust ) {
	    hprintf(wf, "%s\t%d\t%d,%d",hdrs[c].c_str(), colSums[c], k2+1, k1+1);
	    k1 = k2 = nClust + 1;	    
	    break;
	  }
	}
      }
      if ( k1 == nClust ) error("Cannot recognize iBest = %d, nClust = %d", iBest, nClust);
    }
    for(int32_t k=0; k < nClust + nPair; ++k)
      hprintf(wf, "\t%.5lg",Zs[c*(nClust+nPair) + k]/sumZ);
    hprintf(wf, "\n");
  }
  hts_close(wf);    
  
  // free up the memories
  for(int32_t i=0; i < nRow; ++i) {
    free(R[i]);
  }
  //delete[] llks;
  delete[] p0;
  free(pis);
  free(Zs);
  free(Ps);
  free(Ls);
  free(colSums);
  
  return 0;
}

