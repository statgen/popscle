#ifndef __DISCRETE_LOG_HELPER_H
#define __DISCRETE_LOG_HELPER_H

#include "Error.h"

class discrete_log_helper {
 public:
  double offset;
  std::vector<double> values;
  int32_t sz;
  discrete_log_helper(double _offset) : offset(_offset), sz(0) {}
  inline double get(int32_t v) {
    if ( v < 0 )
      error("%lg get(%d) called", offset, v);
    if ( sz <= v ) {
      //fprintf(stderr,"offset = %lg, resizing from %d to %d\n", offset, sz, v+1);
      values.resize((size_t)v+1);
      for(int32_t i = sz; i <= v; ++i)
	values[i] = log(offset + (double)i);
      sz = v+1;
    }
    return values[v];
  }
};

#endif
