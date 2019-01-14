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

#include "sam_ordered_writer.h"

SAMOrderedWriter::SAMOrderedWriter(std::string output_sam_file_name, int32_t _window)
{
    file_name = output_sam_file_name;
    window = _window;
    file = NULL;

    kstring_t mode = {0,0,0};
    kputc('w', &mode);

    if (file_name=="-")
    {
        //do nothing
    }
    else
    {
        if (str_ends_with(file_name, ".sam"))
        {
            //do nothing
        }
        else if (str_ends_with(file_name, ".bam"))
        {
            kputc('b', &mode);
        }
        else if (str_ends_with(file_name, ".ubam"))
        {
            kputs("bu", &mode);
        }	
        else if (str_ends_with(file_name, ".cram"))
        {
            kputs("c", &mode);
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Not a SAM/BAM file: %s\n", __FILE__,__LINE__,__FUNCTION__, file_name.c_str());
            exit(1);
        }
    }

    file = sam_open(file_name.c_str(), mode.s);
    if (file==NULL)
    {
        fprintf(stderr, "[%s:%d %s] Cannot open SAM/BAM file for writing: %s\n", __FILE__,__LINE__,__FUNCTION__, file_name.c_str());
        exit(1);
    }

    hdr = bam_hdr_init();
    linked_hdr = false;
}

/**
 * Duplicates a hdr and sets it.
 */
void SAMOrderedWriter::set_hdr(bam_hdr_t *hdr)
{
    if (this->hdr)
    {
        bam_hdr_destroy(this->hdr);
    }
    this->hdr = bam_hdr_dup(hdr);
    linked_hdr = false;
}

/**
 * Links a header.  This is useful when the VCF file being read has an incomplete header.
 * As the VCF records are read, the incomplete header will be fixed with string type assumptions
 * and the VCF records can be written out without any failure.  The header in the VCF file being
 * written will be incomplete nonetheless and the user should use an alt header when reading the
 * file to bypass the problem.
 */
void SAMOrderedWriter::link_hdr(bam_hdr_t *hdr)
{
    this->hdr = hdr;
    linked_hdr = true;
}

/**
 * Reads next record, hides the random access of different regions from the user.
 */
void SAMOrderedWriter::write_hdr()
{
  if ( sam_hdr_write(file, hdr) != 0 ) {
    error("[E:%s:%d %s] Cannot write the header file", __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }
}

/**
 * Reads next record, hides the random access of different regions from the user.
 */
void SAMOrderedWriter::write(bam1_t *v)
{
    //place into appropriate position in the buffer
    if (window)
    {
        if (!buffer.empty())
        {
            if (bam_get_tid(v)==bam_get_tid(buffer.back()))
            {
                std::list<bam1_t*>::iterator i = buffer.begin();

                for (i=buffer.begin(); i!=buffer.end(); ++i)
                {
                    //equal sign ensures records are kept in original order
                    if (bam_get_pos1(v)>=bam_get_pos1(*i))
                    {
                        buffer.insert(i,v);
                        flush(false);
                        return;
                    }
                }

                if (i==buffer.end())
                {
                    int32_t cutoff_pos1 =  std::max(bam_get_pos1(buffer.front())-window,1);
                    if (bam_get_pos1(v)<cutoff_pos1)
                    {
		      notice("[%s:%d %s] Might not be sorted for window size %d at current record %s:%d < %d (%d [last record] - %d), please increase window size to at least %d.\n", __FILE__,__LINE__,__FUNCTION__, window, bam_get_chrom(hdr, v), bam_get_pos1(v), cutoff_pos1, bam_get_pos1(buffer.front()), window, bam_get_pos1(buffer.front())-bam_get_pos1(v)+1);
                    }
                }

                buffer.insert(i,v);
                flush(false);
            }
            else
            {
                flush(true);
                buffer.push_front(v);
            }
        }
        else
        {
            buffer.push_front(v);
        }

        v = NULL;
    }
    else
    {
      if ( sam_write1(file, hdr, v) < 0 ) {
	error("[E:%s:%d %s] Cannot write the header file", __FILE__, __LINE__, __PRETTY_FUNCTION__);
      }      
    }
}

/**
 * Flush writable records from buffer.
 */
void SAMOrderedWriter::flush()
{
    flush(true);
}

/**
 * Returns record to pool.
 */
void SAMOrderedWriter::store_bam1_into_pool(bam1_t* v)
{
  //bcf_clear(v);
  pool.push_back(v);
  v = NULL;
}

/**
 * Gets record from pool, creates a new record if necessary
 */
bam1_t* SAMOrderedWriter::get_bam1_from_pool()
{
    if(!pool.empty())
    {
        bam1_t* v = pool.front();
        pool.pop_front();
        return v;
    }
    else
    {
        bam1_t* v= bam_init1();
        //bcf_clear(v);
        return v;
    }
};

/**
 * Flush writable records from buffer.
 */
void SAMOrderedWriter::flush(bool force) {
  if (force) {
    while (!buffer.empty()) {
      if ( sam_write1(file, hdr, buffer.back()) < 0 ) {
	  error("[E:%s:%d %s] Cannot write the record", __FILE__, __LINE__, __PRETTY_FUNCTION__);
      }            
      store_bam1_into_pool(buffer.back());
      buffer.pop_back();
    }
  }
  else {
    if (buffer.size()>1) {
      int32_t cutoff_pos1 =  std::max(bam_get_pos1(buffer.front())-window,1);
      
      while (buffer.size()>1) {
	if (bam_get_pos1(buffer.back())<=cutoff_pos1) {
	  if ( sam_write1(file, hdr, buffer.back()) < 0 ) {
	    error("[E:%s:%d %s] Cannot write the record", __FILE__, __LINE__, __PRETTY_FUNCTION__);	    
	  }
	  store_bam1_into_pool(buffer.back());
	  buffer.pop_back();
	}
	else {
	  return;
	}
      }
    }
  }
}

/**
 * Closes the file.
 */
void SAMOrderedWriter::close()
{
    flush(true);
    sam_close(file);
    if (!linked_hdr && hdr) bam_hdr_destroy(hdr);
}
