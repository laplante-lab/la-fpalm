#ifndef KERNEL_CALL_H
#define KERNEL_CALL_H

//#include "kernel.h"

//void CUDAERROR(const char *instr, int lineNumber);

//void cudasafe(cudaError_t err, char* str, int lineNumber);

//*******************************************************************************************

//extern "C" void kernel_fit_one_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float *d_varim, float *d_gainim, float gsig, int subsz, int iter_num, float *d_x, float *d_y, float *d_I, float *d_bg, float *d_C, float *d_ll, int Nfits);

//extern "C" void kernel_fit_two_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float *d_varim, float *d_gainim, float gsig, int subsz, int iter_num, float *d_x, float *d_y, float *d_I, float *d_bg, float *d_s, float *d_C, float *d_ll, int Nfits);

//extern void kernel_MLEFit_sigmaxy_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float *d_varim, float *d_gainim, float PSFSigma, const int sz, const int iterations, float *d_C, float *d_ll, const int Nfits, float *d_x, float *d_y, float *d_I, float *d_bg, float *d_s, float *d_s2);


//void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const	mxArray	*prhs[]) {
void SRsCMOS_MLE( int num_dataset, int subtot, float *subregion, float *subvarim, float *subgainim, int fit_type, float gsig, int iter_num, float *sub_x, float *sub_y, float *sub_photon, float *sub_background, float *CRLB2, float *LL2 );

#endif
