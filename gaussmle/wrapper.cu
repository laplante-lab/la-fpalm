/*
 * \File : wrapper.cu
 * \Author : Peiyi Zhang
 * \Date : November 3rd, 2016
 * \This is a modified version of Fang Huang's wrapper.cu in sCMOS software. 
 * \Copyright @ Purdue University
 */
#include "definitions.h"
#include<device_launch_parameters.h>
#include<device_functions.h>

cudaStream_t stream2;

__global__ void kernel_fit_one(float *d_data, float *d_varim, float *d_gainim, float PSFsigma, int sz, int iterations,
	float *d_x, float *d_y, float *d_I, float *d_bg, float *d_C, float *d_ll, int Nfits);

extern "C" void kernel_fit_one_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float *d_varim, float *d_gainim, float PSFSigma,  int subsz, int iterations, 
        float *d_x, float *d_y, float *d_I, float *d_bg, float *d_C, float *d_ll,int Nfits) 
{
	kernel_fit_one<<<dimGrid, dimBlock, 0, stream2>>>(d_data, d_varim, d_gainim, PSFSigma, subsz, iterations, d_x, d_y, d_I, d_bg, d_C, d_ll, Nfits);
}

__global__ void kernel_fit_two(float *d_data, float *d_varim, float *d_gainim, float PSFsigma, int sz, int iterations,
	float *d_x, float *d_y, float *d_I, float *d_bg, float *d_s, float *d_C, float *d_ll, int Nfits);

extern "C" void kernel_fit_two_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float *d_varim, float *d_gainim, float PSFSigma,  int subsz, int iterations, 
        float *d_x, float *d_y, float *d_I, float *d_bg, float *d_s, float *d_C, float *d_ll,int Nfits) 
{
	kernel_fit_two<<<dimGrid, dimBlock, 0, stream2>>>(d_data, d_varim, d_gainim, PSFSigma, subsz, iterations, d_x, d_y, d_I, d_bg, d_s, d_C, d_ll, Nfits);
}

__global__ void kernel_MLEFit_sigmaxy(float *d_data, float *d_varim, float *d_gainim, float PSFSigma,
	const int sz, const int iterations, float *d_C, float *d_ll, const int Nfits, float *d_x, float *d_y,
	float *d_I, float *d_bg, float *d_s, float *d_s2);

extern void kernel_MLEFit_sigmaxy_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float *d_varim, float *d_gainim, float PSFSigma,const int sz, const int iterations, float *d_C, float *d_ll, const int Nfits, float *d_x, float *d_y,float *d_I, float *d_bg, float *d_s, float *d_s2)
{
	kernel_MLEFit_sigmaxy << <dimGrid, dimBlock, 0, stream2 >> >(d_data, d_varim, d_gainim, PSFSigma, sz, iterations, d_C, d_ll, Nfits, d_x, d_y, d_I, d_bg, d_s, d_s2);
}

