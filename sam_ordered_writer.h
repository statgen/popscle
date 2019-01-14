/* The MIT License

   Copyright (c) 2016 Hyun Min Kang

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

#ifndef SAM_ORDERED_WRITER_H
#define SAM_ORDERED_WRITER_H

#include <cstdlib>
#include <stdint.h>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include "htslib/sam.h"
#include "hts_utils.h"
#include "Error.h"

/**
 * A class for writing ordered SAM/BAM/CRAM files.
 *
 * In some cases, the processing of a file changes the coordinates
 * slightly leading to slight local disorder within a contig,
 * instead of sorting the VCF wholesale, this class buffers the output
 * and sorts locally in a 10K base pair region before writing the records
 * out.
 */
class SAMOrderedWriter
{
    public:

    ///////
    //i/o//
    ///////
    std::string file_name;
    samFile *file;
    bam_hdr_t *hdr;
    bool linked_hdr;

    //buffer for containing records to be written out
    std::list<bam1_t*> buffer; //most recent records in the front
    std::list<bam1_t*> pool;

    int32_t window;

    /**
     * Initialize output file.
     * @output_vcf_file
     * @window - the window to keep variants in buffer to check for local disorder, 0 for no buffering
     * @recycle - keep unused records in a pool for reuse
     */
    SAMOrderedWriter(std::string output_sam_file_name, int32_t _window=0);

    /**
     * Duplicates a hdr and sets it.
     */
    void set_hdr(bam_hdr_t *hdr);

    /**
     * Links a header.  This is useful when the read in VCF header is incomplete.
     */
    void link_hdr(bam_hdr_t *hdr);

    /**
     * Reads next record, hides the random access of different regions from the user.
     */
    void write_hdr();

    /**
     * Reads next record, hides the random access of different regions from the user.
     */
    void write(bam1_t *v);

    /**
     * Gets record from pool, creates a new record if necessary.
     * This is exposed so that the programmer may reuse bam1_t
     * from this class and return to it when writing which is
     * essentially stowing it away in a buffer.
     */
    bam1_t* get_bam1_from_pool();

    /**
     * Returns record to pool.
     */
    void store_bam1_into_pool(bam1_t* v);

    /**
     * Flush writable records from buffer.
     */
    void flush();

    /**
     * Closes the file.
     */
    void close();

    private:

    /**
     * Flush writable records from buffer.
     */
    void flush(bool force);
};

#endif
