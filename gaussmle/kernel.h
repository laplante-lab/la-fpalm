/*
 * \File : mexFunction.cpp
 * \Author : Peiyi Zhang
 * \Date : November 3rd, 2016
 * \This is a modified version of Fang Huang's kernel.h in sCMOS software. 
 * \Copyright @ Purdue University
 */
#include "definitions.h"
#ifndef KERNEL_H
#define KERNEL_H

__global__ void kernel_fit_one(float *d_data, float *d_varim, float *d_gainim, float PSFsigma, int sz, int iterations,
	float *d_x, float *d_y, float *d_I, float *d_bg, float *d_C, float *d_ll, int Nfits);

__global__ void kernel_fit_two(float *d_data, float *d_varim, float *d_gainim, float PSFsigma, int sz, int iterations,
	float *d_x, float *d_y, float *d_I, float *d_bg, float *d_s, float *d_C, float *d_ll, int Nfits);

__global__ void kernel_MLEFit_sigmaxy(float *d_data, float *d_varim, float *d_gainim, float PSFSigma, 
	const int sz, const int iterations, float *d_C, float *d_ll, const int Nfits, float *d_x, float *d_y, 
    float *d_I, float *d_bg, float *d_s, float *d_s2);

// inverse matrix with back substitute method
__device__ void kernel_MatInv(float * m, float * inv, float * crlb, int sz);

#endif
