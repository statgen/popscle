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

#include "utils.h"

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char *delims, std::string& str, uint32_t limit, bool clear, bool collapse)
{
    std::map<char, int32_t> delim_set;

    for (uint32_t i=0; i<strlen(delims); ++i)
    {
        delim_set[delims[i]] = 1;
    }

    if (clear)
    {
        vec.clear();
    }
    const char* tempStr = str.c_str();
    int32_t i=0, lastIndex = str.size()-1;
    std::stringstream token;

    if (lastIndex<0) return;

    uint32_t noTokens = 0;
    bool isDelim = false;
    while (i<=lastIndex)
    {
        isDelim = (delim_set.find(tempStr[i])!=delim_set.end());

        if (!isDelim || noTokens>=limit-1)
        {
            token << tempStr[i];
        }

        if ((isDelim && noTokens<limit-1) || i==lastIndex)
        {
            if (collapse && token.str()!="")
            {
                vec.push_back(token.str());
                ++noTokens;
                token.str("");
            }
        }

        ++i;
    }
};

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char *delims, const char* str, uint32_t limit, bool clear, bool collapse)
{
    std::map<char, int32_t> delim_set;

    for (uint32_t i=0; i<strlen(delims); ++i)
    {
        delim_set[delims[i]] = 1;
    }

    if (clear)
    {
        vec.clear();
    }
    const char* tempStr = str;
    int32_t i=0, lastIndex = strlen(str)-1;
    std::stringstream token;

    if (lastIndex<0) return;

    uint32_t noTokens = 0;
    bool isDelim = false;
    while (i<=lastIndex)
    {
        isDelim = (delim_set.find(tempStr[i])!=delim_set.end());

        if (!isDelim || noTokens>=limit-1)
        {
            token << tempStr[i];
        }

        if ((isDelim && noTokens<limit-1) || i==lastIndex)
        {
            if (collapse && token.str()!="")
            {
                vec.push_back(token.str());
                ++noTokens;
                token.str("");
            }
        }

        ++i;
    }
};

/**
 * Casts a string into int32.  Returns true if successful.
 */
bool str2int32(std::string& s, int32_t& i)
{
    const char* start = s.c_str();
    char *end = 0;
    i = std::strtol(s.c_str(), &end, 10);
    return (end!=start);
};

/**
 * Casts a string into uint32.  Returns true if successful.
 */
bool str2uint32(std::string& s, uint32_t& i)
{
    const char* start = s.c_str();
    char *end = 0;
    i = std::strtoul(s.c_str(), &end, 10);
    return (end!=start);
};

/**
 * Appends cuurent working directoy to a path.
 * Returns true if successful.
 */
bool append_cwd(std::string& path)
{
    if (path.size()>0 && path.c_str()[0]!='/')
    {
        char cwd[1024];
        if (getcwd(cwd, sizeof(cwd))!=NULL)
        {
            std::string cwd_path(cwd);
            path = cwd_path + "/" + path;

            return true;
        }
    }

    return false;
};


unsigned int str_hash(const char* s, unsigned int seed)
{
  unsigned int hash = seed;
  while (*s) {
    hash = hash * 101  +  *s++;
  }
  return hash;
}

