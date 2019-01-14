/* The MIT License

   Copyright (c) 2015 Adrian Tan <atks@umich.edu>

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

#ifndef REFERENCE_SEQUENCE_H
#define REFERENCE_SEQUENCE_H

#include "utils.h"
#include "hts_utils.h"

/**
 * A Reference Sequence object wrapping htslib's faidx.
 * This allows for buffered reading of seqeunces.
 */
class ReferenceSequence
{
    public:

    //reference file and index
    std::string ref_fasta_file;
    faidx_t *fai;

    uint32_t buffer_size;
    uint32_t buffer_size_mask;
    uint32_t window_size;

    //main buffer sequence
    std::vector<char> seq;

    //secondary buffer sequence for wrap around cases of the sequence.
    std::vector<char> cseq;

    int32_t tid;
    uint32_t beg0, end0; // index of P for the start and end position in 0 base coordinates

    std::string chrom;
    uint32_t gbeg1;

    //  beg0 points to the beginning of the buffered sequence
    //  end0 points to the position after the end of the buffered sequence
    //
    //      end0
    //       |
    //       V
    //  +--+--+--+--+--+--+--+--+--+--+--+--+
    //  |  |  |  |  |  |  |  |  |  |  |  |  |
    //  +--+--+--+--+--+--+--+--+--+--+--+--+
    //       ^
    //       |
    //      beg0
    //
    //  when beg0==end0, buffered sequence is empty
    //  gbeg1 - is the coordinate of the genome position that beg0 represents.
    //
    //  invariance:  buffered sequence is always continuous

    int32_t debug;

    public:

    /**
     * Constructor.
     *
     * @k - size of buffered sequence is 2^k.
     */
    ReferenceSequence(std::string& ref_fasta_file, uint32_t k=10, uint32_t window_size=256);

    /**
     * Fetches the number of sequences.
     */
    int32_t fetch_nseq();

    /**
     * Fetch name of the ith sequence.
     */
    std::string fetch_iseq_name(int32_t i);

    /**
     * Fetch length of sequence seq.
     */
    int32_t fetch_seq_len(std::string& seq);

    int32_t fetch_seq_len(const char* chrom) {
      std::string s(chrom);
      return fetch_seq_len(s);
    }

    

    /**
     * Get a base.
     */
    char fetch_base(const char* chrom, int32_t pos1);

    /**
     * Get a base.
     */
    char fetch_base(std::string& chrom, int32_t pos1);

    /**
     * Fetches sequence chrom:beg1-end1.
     *
     * Retrieved sequence is in seq with the length of n.
     */
    //void fetch_seq(std::string& chrom, int32_t start1, int32_t end1, char* seq, int32_t n);

    /**
     * Fetches sequence chrom:beg1-end1.
     */
    void fetch_seq(std::string& chrom, int32_t beg1, int32_t end1, std::string& seq);

    /**
     * Fetches sequence chrom:beg1-end1.
     */
    void fetch_seq(const char* chrom, int32_t beg1, int32_t end1, std::string& seq);

    /**
     * Fetches sequence chrom:beg1-end1.
     */
    char* fetch_seq(const char* chrom, int32_t beg1, int32_t end1);

    /**
     * Fetches sequence chrom:beg1-end1.
     */
    char* fetch_seq(char* chrom, int32_t beg1, int32_t end1);

    /**
     * Fetches sequence chrom:beg1-end1.
     */
    char* fetch_seq(std::string& chrom, int32_t beg1, int32_t end1);

    private:

    /**
     * Overloads subscript operator for accessing buffered sequence positions.
     */
    char& operator[] (const int32_t i);

    /**
     * Returns the maximum size of the buffered sequence.
     */
    uint32_t max_size();

    /**
     * Returns the size of the buffered sequence.
     */
    uint32_t size();

    /**
     * Checks if buffer is empty
     */
    bool is_empty();

    /**
     * Set reference fasta file.
     */
    void set_reference(std::string& ref_fasta_file);

    /**
     * Set debug.
     */
    void set_debug(int32_t debug);

    /**
     * Sets tid.
     */
    void set_tid(uint32_t tid);

    /**
     * Gets tid.
     */
    uint32_t get_tid();

    /**
     * Sets chrom.
     */
    void set_chrom(std::string& chrom);

    /**
     * Gets chrom.
     */
    std::string get_chrom();

    /**
     * Gets window_size.
     */
    uint32_t get_window_size();

    /**
     * Converts gpos1 to index in P.
     */
    uint32_t g2i(uint32_t gpos1);

    /**
     * Returns the difference between 2 buffer positions
     */
    uint32_t diff(uint32_t i, uint32_t j);

    /**
     * Print buffered sequence state.
     */
    void print_state();

    /**
     * Print buffered sequence state extent.
     */
    void print_state_extent();
};

#endif
