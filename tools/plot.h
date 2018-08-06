#ifndef PLOT_H
#define PLOT_H

#include "../SRscmos.h"

struct fe_data{
  const char* output_file_path;
  float* data;
  int size;
};

struct im_data{
  const char* output_file_path;
  float* xtot;
  float* ytot;
  float* bgtot;
  float* photot;
  float* crlbytot;
  float* lltot;
  int size;
};

struct un_data{
  const char* output_file_path;
  float* crlbxtot;
  float* crlbytot;
  int size;
};

void plot( const char* output_path, SRscmos *sr ); 

#endif
