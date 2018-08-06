#ifndef SRHIST_COLOR_H
#define SRHIST_COLOR_H

double* imstretch_linear(int image_size, double* im, int low_in, int high_in, int low_out, int high_out);

void srhist_color(int sz, int zm, int size, float* xtot, float* ytot, float* ttot, int segnum, double* rch, double* gch, double* bch );
#endif
