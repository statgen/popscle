#ifndef BCF_CHUNKED_READER_H
#define BCF_CHUNKED_READER_H

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <list>
#include <map>
#include <queue>
#include "hts_utils.h"
#include "genome_interval.h"
#include "interval_tree.h"
#include "Error.h"
#include "genomeChunk.h"
#include "genomeLoci.h"

extern "C" {
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"
}

/**
 * A class for reading chunked VCF/BCF files.
 *
 */

class BCFChunkedReader
{
public:
  
  ///////
  //i/o//
  ///////
  
  // information needed if pattern is used
  genomeChunk chunk;
  vcfFile *file;
  bcf_hdr_t *hdr;
  hts_idx_t *idx;
  tbx_t *tbx;
  hts_itr_t *itr;
  bcf1_t *v;
  //bcf1_t *buffer_v;
  
  //for control
  htsFormat ftype;
  //bool intervals_present;
  //bool index_loaded;
  //bool random_access_enabled;
  
  //list of intervals
  genomeLoci target_intervals;
  
  //for storing unused bcf records
  //std::list<bcf1_t*> pool;
  
  //shared objects for string manipulation
  kstring_t s;
  
  /**
   * Initialize files and intervals.
   *
   * @input_vcf_file_name     name of the input VCF file
   * @intervals          list of intervals, if empty, all records are selected.
   */
  void init(genomeLoci* pIntervals);
  
  void init(const char* input_pattern_name, const char* input_ref_file, const char* input_interval_list, int32_t unit, genomeLoci* pIntervals);
  
  bool open_current_file();
  
  BCFChunkedReader() : file(NULL), hdr(NULL), idx(NULL), tbx(NULL), itr(NULL), v(NULL) {}
  
  BCFChunkedReader(const char* input_vcf_file_name, genomeLoci* pIntervals = NULL);
  
  BCFChunkedReader(const char* input_pattern_name, const char* input_ref_file, int32_t unit_chunk, genomeLoci* pIntervals = NULL);
  
  BCFChunkedReader(const char* input_pattern_name, const char* input_interval_list, genomeLoci* pIntervals = NULL);
  
  BCFChunkedReader(const char* input_pattern_name, const char* input_ref_file, const char* input_interval_list, int32_t unit, genomeLoci* pIntervals);
  
  
  /**
   * Returns next vcf record.
   */
  bool read(bcf1_t *v);
  
  /**
   * Initialize next interval.
   * Returns false only if all intervals are accessed.
   */
  bool initialize_current_interval();
  
  //    /**
  //     * Returns next set of vcf records at a start position.
  //     * Note that this function should never be used in conjunction with read(bcf1_t *v)
  //     */
  //    bool read_next_position(std::vector<bcf1_t *>& vs);
  
  /**
   * Gets sequence name of a record.
   */
  const char* get_seqname(bcf1_t *v);
  
  /**
   * Checks if index is loaded.
   */
  //bool is_index_loaded();
  
  /**
   * Gets bcf header.
   */
  bcf_hdr_t* get_hdr();
  
  /**
   * Closes the file.
   */
  void close(bool destroy_hdr = true);
  
  bool load_index();
  
  bool jump_to(const char* chr = NULL, int32_t pos1 = INT_MAX);
};

#endif
