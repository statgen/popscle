//#include <iostream>
//#include <string>
//#include <map>
//#include <vector>
//#include <ctime>
//#include <cstdlib>
//#ifdef _OPENMP
//  #include <omp.h>
//#else
//  #define omp_get_thread_num() 0
//#endif  

#include "cramore.h"
#include "commands.h"

//#include "htslib/sam.h"
//#include "htslib/faidx.h"
//#include "htslib/kseq.h"
//#include "htslib/khash.h"
//#include "dropseq.h"
#include "utils.h"

KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;


int32_t cmdCramDemuxlet(int32_t argc, char** argv);
int32_t cmdCramFreemuxlet(int32_t argc, char** argv);
int32_t cmdCramMuxPileup(int32_t argc, char** argv);
int32_t cmdCramDigitalPileup(int32_t argc, char** argv);
int32_t cmdCramSimuxlet(int32_t argc, char** argv);



int32_t cmdScMultinomEM(int32_t argc, char** argv);
int32_t cmdScMultinomGibbs(int32_t argc, char** argv);
int32_t cmdScMapSTAMPs(int32_t argc, char** argv);
int32_t cmdScKallistoCount(int32_t argc, char** argv);

int32_t cmdGtfUtil(int32_t argc, char** argv);
int32_t cmdPlpMakeDGEMatrix(int32_t argc, char** argv);
int32_t cmdPlpFindMuxClust(int32_t argc, char** argv);


int32_t main(int32_t argc, char** argv) {
  commandList cl;

  BEGIN_LONG_COMMANDS(longCommandlines)
    
    LONG_COMMAND_GROUP("Single cell analysis", NULL)
    LONG_COMMAND("demuxlet", &cmdCramDemuxlet, "Deconvolute sample identify of droplet-based sc-RNAseq")
    LONG_COMMAND("freemuxlet", &cmdCramFreemuxlet, "Genotype-free deconvolution of sc-RNAseq")
    //LONG_COMMAND("mux-pileup", &cmdCramMuxPileup, "Produce pileup of dsc-RNAseq")
    LONG_COMMAND("dsc-pileup", &cmdCramDigitalPileup, "Produce pileup of dsc-RNAseq")    
    //LONG_COMMAND("simuxlet",   &cmdCramSimuxlet,  "Simulate multiplexed dsc-RNAseq droplets")
    //LONG_COMMAND("kallisto-count", &cmdScKallistoCount, "Produce digital expression matrix from kallisto-aligned sequence reads")
    //LONG_COMMAND("sc-map-stamps", &cmdScMapSTAMPs, "Produce STAMP-map from DropSeq FASTQ files")
    //LONG_COMMAND("sc-multinom-em", &cmdScMultinomEM, "Multinomial EM clustering of single cell types")
    //LONG_COMMAND("sc-multinom-gibbs", &cmdScMultinomGibbs, "Multinomial Gibbs sampling of single cell types")
    //LONG_COMMAND("plp-find-mux-clust", &cmdPlpFindMuxClust, "(Experimental) Perform tSNE clustering on mux-seq data")

    //LONG_COMMAND_GROUP("Other tools", NULL)
    //LONG_COMMAND("gtf-util", &cmdGtfUtil, "A simple GTF utility")
    //LONG_COMMAND("plp-make-dge-matrix", &cmdPlpMakeDGEMatrix, "Make Digital Expression Matrix from Digital Pileups")
  END_LONG_COMMANDS();

  cl.Add(new longCommands("Available Commands", longCommandlines));
  
  if ( argc < 2 ) {
    printf("[freemuxlet] -- Genotyping-free deconvolution tool for multiplexed single cell experiment over multiple individuals\n\n");
    fprintf(stderr, " Copyright (c) 2009-2019 by Hyun Min Kang and Fan Zhang\n");
    fprintf(stderr, " Licensed under the Apache License v2.0 http://www.apache.org/licenses/\n\n");    
    fprintf(stderr, "To run a specific command      : %s [command] [options]\n",argv[0]);
    fprintf(stderr, "For detailed instructions, run : %s --help\n",argv[0]);        
    cl.Status();
    return 1;
  }
  else {
    if ( strcmp(argv[1],"--help") == 0 ) {
      cl.HelpMessage();
    }
    else {
      return cl.Read(argc, argv);
    }
  }
  return 0;
}
