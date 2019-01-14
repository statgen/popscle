#ifndef __CRAMORE_H
#define __CRAMORE_H

#include <getopt.h>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <set>
#include <ctime>
#include <cmath>

// Hyun's codes
#include "params.h"
#include "Error.h"
#include "PhredHelper.h"

// Adrian's codes
//#include "genome_interval.h"
#include "hts_utils.h"

extern "C" {
  size_t hts_realloc_or_die(unsigned long, unsigned long, unsigned long, unsigned long, int, void**, char const*);
}

// bcftools's code
//#include "filter.h"
//#include "genomeChunk.h"
//#include "genomeLoci.h"

// htslib's code

#endif
