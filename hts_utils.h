/* The MIT License

   Copyright (c) 2013 Adrian Tan <atks@umich.edu>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#ifndef HTS_UTILS_H
#define HTS_UTILS_H

#include <string>
#include <vector>
#include <cstdlib>
//#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <map>
#include <queue>
#include <sys/stat.h>

extern "C" {
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"
}

#include "utils.h"
#include "genome_interval.h"

/**********
 *FAI UTILS
 **********/
typedef struct {
    int32_t line_len, line_blen;
    int64_t len;
    uint64_t offset;
} faidx1_t;

KHASH_MAP_INIT_STR(s, faidx1_t)

struct __faidx_t {
    BGZF *bgzf;
    int n, m;
    char **name;
    khash_t(s) *hash;
};

/**
 * An alternate sequence fetcher for upper case sequence.
 */
char *faidx_fetch_uc_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len);

/**********
 *HTS UTILS
 **********/

/**
 * Checks file extension for use in writing files.
 */
bool str_ends_with(std::string& file_name, const char* ext);

/**************
 *BAM HDR UTILS
 **************/

/**
 * Copies contigs found in bam header to bcf header.
 */
void bam_hdr_transfer_contigs_to_bcf_hdr(const bam_hdr_t *sh, bcf_hdr_t *vh);

/**
 * Get number of sequences.
 */
#define bam_hdr_get_n_targets(h) ((h)->n_targets)

/**
 * Get list of sequence names.
 */
#define bam_hdr_get_target_name(h) ((h)->target_name)

/**
 * Get list of sequence lenghts.
 */
#define bam_hdr_get_target_len(h) ((h)->target_len)

/**********
 *BAM UTILS
 **********/

//used when a base on a read is not aligned to the human genome reference
//in particular for soft clips and insertions
#define BAM_READ_INDEX_NA -1

/**
 * Gets the chromosome name of the tid.
 */
#define bam_get_chrom(h, s) ((h)->target_name[(s)->core.tid])

#define bam_get_chromi(h, i) ((h)->target_name[i])

/**
 * Gets the 1 based start position of the first mapped base in the read.
 */
#define bam_get_pos1(s) ((s)->core.pos + 1)

/**
 * Gets the 0 based start position of the first mapped base in the read.
 */
#define bam_get_pos0(s) ((s)->core.pos)

/**
 * Gets the end position of the last mapped base in the read.
 */
int32_t bam_get_end_pos1(bam1_t *srec);

/**
 * Gets the template ID of the read.
 */
#define bam_get_tid(s) ((s)->core.tid)

/**
 * Gets the template ID of the paired read.
 */
#define bam_get_mtid(s) ((s)->core.mtid)

/**
 * Gets the start position of the first mapped base in the read.
 */
#define bam_get_mpos1(s) ((s)->core.mpos)

/**
 * Gets the read sequence from a bam record
 */
void bam_get_seq_string(bam1_t *s, kstring_t *seq);

/**
 * Gets the base qualities from a bam record.
 */
void bam_get_qual_string(bam1_t *s, kstring_t *qual);

/**
 * Gets the cigar sequence from a bam record.
 */
#define bam_get_n_cigar_op(b) ((b)->core.n_cigar)

/**
 * Gets the cigar from a BAM record.
 */
void bam_get_cigar_string(bam1_t *s, kstring_t *cigar);

/**
 * Gets the cigar string from a bam record.
 */
void bam_get_cigar_expanded_string(bam1_t *s, kstring_t *cigar_string);

/**
 * Is this sequence the first read?
 */
#define bam_is_fread1(s) (((s)->core.flag&BAM_FREAD1) != 0)

/**
 * Is this sequence the second read?
 */
#define bam_is_fread2(s) (((s)->core.flag&BAM_FREAD2) != 0)

/**
 * Gets the base in the read that is mapped to a genomic position.
 * Returns -1 if it does not exists.
 */
void bam_get_base_and_qual(bam1_t *srec, uint32_t pos, char& base, char& qual, int32_t& rpos);

/**
 * Gets the base in the read that is mapped to a genomic position.
 * Returns true if successful, false if the location is a deletion or not on the read.
 */
void bam_get_base_and_qual_and_read_and_qual(bam1_t *s, uint32_t pos, char& base, char& qual, int32_t& rpos, kstring_t* readseq, kstring_t* readqual);

/**
 * Converts base encoding to character.
 */
#define bam_base2char(b) ("XACXGXXXTXXXXXXN"[(b)])

/**
 * Gets flag.
 */
#define bam_get_flag(s) ((s)->core.flag)

/**
 * Get map quality.
 */
#define bam_get_mapq(s) ((s)->core.qual)

/**
 * Get tid - e.g. chromosome id.
 */
#define bam_get_tid(s) ((s)->core.tid)

/**
 * Get read length.
 */
#define bam_get_l_qseq(s) ((s)->core.l_qseq)

/**
 * Prints a bam.
 */
void bam_print(bam_hdr_t *h, bam1_t *s);

/**************
 *BCF HDR UTILS
 **************/

/**
 * Copies contigs found in bcf header to another bcf header.
 */
void bcf_hdr_transfer_contigs(const bcf_hdr_t *sh, bcf_hdr_t *vh);

/**
 * Extracts sequence length by rid.
 */
int32_t* bcf_hdr_seqlen(const bcf_hdr_t *hdr, int32_t *nseq);

/**
 * Get samples from bcf header
 */
#define bcf_hdr_get_samples(h) ((h)->samples)

/**
 * Get ith sample name from bcf header
 */
#define bcf_hdr_get_sample_name(h, i) ((h)->samples[i])

/**
 * Get number of samples in bcf header
 */
int32_t bcf_hdr_get_n_sample(bcf_hdr_t *h);

/**
 * Reads header of a VCF file and returns the bcf header object.
 * This wraps around vcf_hdr_read from the original htslib to
 * allow for an alternative header file to be read in.
 *
 * this searches for the alternative header saved as <filename>.hdr
 * If the VCF files is BCF, any alternative header is ignored.
 */
bcf_hdr_t *bcf_alt_hdr_read(htsFile *fp);

/**
 * Subsets a record by samples.
 * @imap - indices the subsetted samples
 */
int bcf_hdr_subset_samples(const bcf_hdr_t *h, bcf1_t *v, std::vector<int32_t>& imap);

/**
 * Help function for adding a header with a backup tag name.
 * If the <tag> is already present, a new tag is attempted
 * in the format <tag>_1 to <tag>_9.  If <tag>_9 failed,
 * the function will not add any new tag and will return
 * an empty string.
 *
 * Returns the tag that was inserted or updated.
 */
std::string bcf_hdr_append_info_with_backup_naming(bcf_hdr_t *h, std::string tag, std::string number, std::string type, std::string description, bool rename);

/**********
 *BCF UTILS
 **********/

/**
 * Gets number of expected genotypes from number of allelles for a ploidy of 2.
 */
#define bcf_an2gn(n) (((n+1)*n)>>1)

/**
 * n choose r.
 */
uint32_t choose(uint32_t n, uint32_t r);

/**
 * Gets number of genotypes from number of alleles and ploidy.
 */
uint32_t bcf_ap2g(uint32_t no_allele, uint32_t no_ploidy);

/**
 * Gets number of genotypes from number of alleles and genotypes.
 */
uint32_t bcf_ag2p(uint32_t no_alleles, uint32_t no_genotypes);

/**
 * Gets index from 2 alleles and 2 ploidy.
 */
#define bcf_2g2c(i,j) ((i)+((((j)+1)*(j))>>1))

/**
 * Gets number of genotypes from number of alleles and ploidy.
 */
uint32_t bcf_g2i(std::string genotype);

/**
 * Returns true if a is before b, false otherwise.
 */
bool bcf_is_in_order(bcf1_t *a, bcf1_t *b);

/**
 * Gets a string representation of a variant.
 */
void bcf_variant2string(bcf_hdr_t *h, bcf1_t *v, kstring_t *var);

/**
 * Gets a sorted string representation of a variant.
 */
void bcf_variant2string_sorted(bcf_hdr_t *h, bcf1_t *v, kstring_t *var);

/**
 * Gets a string representation of the alleles of a variant.
 */
void bcf_alleles2string(bcf_hdr_t *h, bcf1_t *v, kstring_t *var);

/**
 * Gets a sorted string representation of the alleles of a variant.
 */
void bcf_alleles2string_sorted(bcf_hdr_t *h, bcf1_t *v, kstring_t *var);

/**
 * Prints a VCF record to STDERR.
 */
void bcf_print(bcf_hdr_t *h, bcf1_t *v);

/**
 * Prints a VCF record in compact string representation to STDERR.
 */
void bcf_print_lite(bcf_hdr_t *h, bcf1_t *v);

/**
 * Prints a VCF record in compact string representation to STDERR with alleles sorted.
 */
void bcf_print_lite_sorted(bcf_hdr_t *h, bcf1_t *v);

/**
 * Prints a VCF record in compact string representation to STDERR with a new line.
 */
void bcf_print_liten(bcf_hdr_t *h, bcf1_t *v);

/**
 * Get number of chromosomes
 */
#define bcf_get_n_chrom(h) ((h)->n[BCF_DT_CTG])

/**
 * Get chromosome name
 */
const char* bcf_get_chrom(bcf_hdr_t *h, bcf1_t *v);

/**
 * Set chromosome name
 */
void bcf_set_chrom(bcf_hdr_t *h, bcf1_t *v, const char* chrom);

/**
 * Get RID
 */
#define bcf_get_rid(v) ((v)->rid)

/**
 * Set RID
 */
#define bcf_set_rid(v, c) ((v)->rid=(c))

/**
 * Check if variant is passed
 */
bool bcf_is_passed(bcf_hdr_t *h, bcf1_t *v);

/**
 * Get 0-based position
 */
#define bcf_get_pos0(v) ((v)->pos)

/**
 * Set 0-based position
 */
#define bcf_set_pos0(v, p) ((v)->pos = (p))

/**
 * Get 1-based position
 */
#define bcf_get_pos1(v) ((v)->pos+1)

/**
 * Get 1-based end position
 */
#define bcf_get_end_pos1(v) ((v)->pos + strlen((v)->d.allele[0]))

/**
 * Set 1-based position
 */
#define bcf_set_pos1(v, p) ((v)->pos = (p)-1)

/**
 * Get id
 */
#define bcf_get_id(v) ((v)->d.id)

/**
 * Set id.
 */
void bcf_set_id(bcf1_t *v, char* id);

/**
 * Get allele
 */
#define bcf_get_allele(v) ((v)->d.allele)

/**
 * Get n_alleles of a bcf record
 */
#define bcf_get_n_allele(v) ((v)->n_allele)

/**
 * Get reference base for a SNP, assumes the record is a biallelic SNP
 */
#define bcf_get_snp_ref(v) ((v)->d.allele[0][0])

/**
 * Get alternative base for a SNP, assumes the record is a biallelic SNP
 */
#define bcf_get_snp_alt(v) ((v)->d.allele[1][0])

/**
 * Get reference allele
 */
#define bcf_get_ref(v) ((v)->d.allele[0])

/**
 * Get alternative allele
 */
#define bcf_get_alt(v, i) ((v)->d.allele[i])

/**
 * Get variant type
 */
#define bcf_get_var_type(v) ((v)->d.var_type)

/**
 * Get qual
 */
#define bcf_get_qual(v) ((v)->qual)

/**
 * Get n_filt of a bcf record
 */
#define bcf_get_n_filter(v) ((v)->d.n_flt)

/**
 * Get ith format name
 */
#define bcf_get_format(h, v, i) (h)->id[BCF_DT_ID][(v->d.fmt)[i].id].key
/**
 * Get number of samples in bcf record
 */
#define bcf_get_n_sample(v) ((v)->n_sample)

const char* bcf_hdr_sample_id(bcf_hdr_t* h, int32_t idx);
int32_t bcf_hdr_sample_index(bcf_hdr_t* h, const char* id);

/**
 * Set number of samples in bcf record
 */
#define bcf_set_n_sample(v, n) ((v)->n_sample = (n))

/**
 * Get qual
 */
#define bcf_get_qual(v) ((v)->qual)

/**
 * Set qual
 */
#define bcf_set_qual(v, q) ((v)->qual = (q))

void hprintf(htsFile* fp, const char * msg, ...);

class GenomeInterval;
void parse_intervals(std::vector<GenomeInterval>& intervals, std::string interval_list, std::string interval_string);

std::string bam_hdr_get_sample_name(bam_hdr_t* hdr);

int32_t bam_get_unclipped_start(bam1_t* b);
int32_t bam_get_unclipped_end(bam1_t* b);
int32_t bam_get_clipped_end(bam1_t* b);

bool same_hrecs(bcf_hdr_t* dst_hdr, bcf_hrec_t* dst, bcf_hdr_t* src_hdr, bcf_hrec_t* src);

char *samfaipath(const char *fn_ref);

bam_hdr_t *sam_hdr_sanitise(bam_hdr_t *h);

//bam_hdr_t* bam_hdr_merge(std::vector<bam_hdr_t*>& hdrs);

#endif
