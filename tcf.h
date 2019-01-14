#ifndef HTSLIB_MCF_H
#define HTSLIB_MCF_H

#include <stdint.h>
#include "hts.h"
#include "kstring.h"
#include "hts_defs.H"
#include "vcf.h"

typedef struct {
  int32_t rid;
  int32_t pos;
  int32_t rlen;
  float qual;
  uint32_t n_info:16, n_allele:16;
  uint32_t n_fmt:8
} mcf1_t;
  
#endif
