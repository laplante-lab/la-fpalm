/*
 * \File : kernel.cu
 * \Author : Peiyi Zhang
 * \Date : November 3rd, 2016
 * \This file contains all of the cuda kernels
 * \The first 2 kernels are modified version of Fang Huang's kernel.cu in sCMOS software. And the third kernel is a modified version of Keith Lidke's GPUgaussMLEv2.cu.
 * \Copyright @ Purdue University
 */
#include "cuda_runtime.h"
#include "kernel.h"
#include "definitions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<device_launch_parameters.h>
//#include<device_functions.h>
#include<math.h>

// inverse matrix with back substitute method
__device__ void kernel_MatInv(float * m, float * inv, float * crlb, int sz) {
	int ii, jj, kk, num, b;
	float tmp1 = 0;
	float yy[25];

	for (jj = 0; jj < sz; jj++) {
		//calculate upper matrix
		for (ii = 0; ii <= jj; ii++)
			//deal with ii-1 in the sum, set sum(kk=0->ii-1) when ii=0 to zero
		if (ii>0) {
			for (kk = 0; kk <= ii - 1; kk++) tmp1 += m[ii + kk*sz] * m[kk + jj*sz];
			m[ii + jj*sz] -= tmp1;
			tmp1 = 0;
		}

		for (ii = jj + 1; ii<sz; ii++)
		if (jj>0) {
			for (kk = 0; kk <= jj - 1; kk++) tmp1 += m[ii + kk*sz] * m[kk + jj*sz];
			m[ii + jj*sz] = (1 / m[jj + jj*sz])*(m[ii + jj*sz] - tmp1);
			tmp1 = 0;
		}
		else { m[ii + jj*sz] = (1 / m[jj + jj*sz])*m[ii + jj*sz]; }
	}

	tmp1 = 0;

	for (num = 0; num < sz; num++) {
		// calculate yy
		if (num == 0) yy[0] = 1;
		else yy[0] = 0;

		for (ii = 1; ii < sz; ii++) {
			if (ii == num) b = 1;
			else b = 0;
			for (jj = 0; jj <= ii - 1; jj++) tmp1 += m[ii + jj*sz] * yy[jj];
			yy[ii] = b - tmp1;
			tmp1 = 0;
		}

		// calculate inv
		inv[sz - 1 + num*sz] = yy[sz - 1] / m[(sz - 1) + (sz - 1)*sz];

		for (ii = sz - 2; ii >= 0; ii--) {
			for (jj = ii + 1; jj < sz; jj++) tmp1 += m[ii + jj*sz] * inv[jj + num*sz];
			inv[ii + num*sz] = (1 / m[ii + ii*sz])*(yy[ii] - tmp1);
			tmp1 = 0;
		}
	}

	if (crlb) for (ii = 0; ii < sz; ii++) crlb[ii] = inv[ii*sz + ii];

	return;

}

__global__ void kernel_fit_one(float *d_data, float *d_varim, float *d_gainim,
	float PSFsigma, int sz, int iterations, float *d_x, float *d_y, float *d_I, float *d_bg, float *d_C, float *d_ll, int Nfits)
{
	//declare variables and set initial value
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int vnum = 4, count = 0, count1 = 0, countx = 0, county = 0, xx = 0, yy = 0;
	float rt = 0.0, PSFx = 0.0, PSFy = 0.0, model = 0.0, data = 0.0, r_data, r_varim, r_gainim, xlim = 1.0, ylim = 1.0, Ilim = 200, bglim = 2.0, deriv1[4] = { 0 }, deriv2[4] = { 0 }, sum = 0.0, tmp = 0.0, I = 0.0, bg = 1e10, x, y, llr = 0;

	__shared__ float s_MM[16], s_CRLB[4], s_invM[16];
	__shared__ float  s_arr[441], s_arrx[441], s_arry[441];//Matrices to store temporal data before reduction

	//initialization of shared memory
	if (threadIdx.x < 16) { s_MM[threadIdx.x] = 0.0f; s_invM[threadIdx.x] = 0.0f; }
	if (threadIdx.x < 4){ s_CRLB[threadIdx.x] = 0.0f; }
	s_arrx[threadIdx.x] = 0.0f;
	s_arry[threadIdx.x] = 0.0f;
	s_arr[threadIdx.x] = 0.0f;

	//initial guess of x and y
	x = (float)sz / 2;
	y = (float)sz / 2;

	//*******************************************************************************************
	//load data, variancemap and gainmap from global memory to register 
	r_data = d_data[idx];
	r_varim = d_varim[idx];
	r_gainim = d_gainim[idx];

	//loop through each pixel of one subregion to refine bg and photon initial guess. scmos noise ignored
	float norm = 1.0f / 2.0f / PSFsigma / PSFsigma;
	for (xx = 0; xx < sz; xx++)
	{
		for (yy = 0; yy < sz; yy++)
		{
			//initialize
			sum = 0;
			tmp = 0;

			//sync threads to make sure threads are in the same loop
			__syncthreads();

			//C++ use Row-major order. threadIdx%sz represents the column index here, which is row index in Matlab of a subregion 
			s_arry[threadIdx.x] = exp(-pow((float)(((int)(threadIdx.x % sz)) - xx), 2)*norm)*exp(-pow((float)(yy - ((int)(threadIdx.x / sz))), 2)*norm);
			s_arrx[threadIdx.x] = s_arry[threadIdx.x] * r_data;

			//sync threads. So that s_arrx and s_arry is updated before sum
			__syncthreads();

			//Calculate sum and tmp. Use parallel reduction algorithm to speed up
			for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
			{
				if (threadIdx.x < s)
				{
					s_arrx[threadIdx.x] += s_arrx[threadIdx.x + s];
					s_arry[threadIdx.x] += s_arry[threadIdx.x + s];
				}
				//sync threads before next loop
				__syncthreads();

				if (threadIdx.x == 0)
				{
					s_arrx[threadIdx.x] += s_arrx[2 * s];//write the sum of elements in s_arrx[] in s_arrx[0]
					s_arry[threadIdx.x] += s_arry[2 * s];//write the sum of elements in s_arry[] in s_arry[0]
				}
			}

			//sync threads to wait for s_arrx[0] and s_arry[0] being calculated before load them to registers of each thread
			__syncthreads();

			tmp = s_arrx[0];
			sum = s_arry[0];

			tmp = tmp / sum;

			if (I < tmp){
				x = (float)xx;
				y = (float)yy;
				I = tmp;
			}

			bg = min(bg, tmp);

		}
	}

	//constraint of background and photon
	bg = max(bg, 2);
	I = max(0.0, (I - bg) * 2 * pi*PSFsigma*PSFsigma);

	//*******************************************************************************************
	// iteration starts
	for (count = 0; count < iterations; count++)
	{
		//initialize
		if (threadIdx.x < 16) { s_MM[threadIdx.x] = 0; }

		PSFx = 1.0f / 2.0f*(erf(float(((int)(threadIdx.x % sz)) - x + 0.5f)*sqrt(norm)) - erf(float(((int)(threadIdx.x % sz)) - x - 0.5f)*sqrt(norm)));
		PSFy = 1.0f / 2.0f*(erf(float(((int)(threadIdx.x / sz)) - y + 0.5f)*sqrt(norm)) - erf(float(((int)(threadIdx.x / sz)) - y - 0.5f)*sqrt(norm)));
		model = bg + I*PSFx*PSFy + r_varim / r_gainim / r_gainim;
		data = r_data + r_varim / r_gainim / r_gainim;

		//calculate derivatives
		deriv1[0] = -I / sqrt(2.0f*pi) / PSFsigma*(exp(-1.0f / 2.0f*pow((float(((int)(threadIdx.x % sz)) + 0.5f - x) / PSFsigma), 2.0f)) - exp(-1.0f / 2.0f*pow((float(((int)(threadIdx.x % sz)) - 0.5f - x) / PSFsigma), 2.0f)))*PSFy;
		deriv1[1] = -I / sqrt(2.0f*pi) / PSFsigma*(exp(-1.0f / 2.0f*pow((float(((int)(threadIdx.x / sz)) + 0.5f - y) / PSFsigma), 2.0f)) - exp(-1.0f / 2.0f*pow((float(((int)(threadIdx.x / sz)) - 0.5f - y) / PSFsigma), 2.0f)))*PSFx;
		deriv1[2] = PSFx * PSFy;
		deriv1[3] = 1.0f;

		deriv2[0] = -I / sqrt(2.0f*pi) / pow(PSFsigma, 3)*(float(((int)(threadIdx.x % sz)) + 0.5f - x)*exp(-1.0f / 2.0f*pow((float(((int)(threadIdx.x % sz)) + 0.5f - x) / PSFsigma), 2.0f)) - float(((int)(threadIdx.x % sz)) - 0.5f - x)*exp(-1.0f / 2.0f*pow(float(((int)(threadIdx.x % sz)) - 0.5f - x) / PSFsigma, 2.0f)))*PSFy;
		deriv2[1] = -I / sqrt(2.0f*pi) / pow(PSFsigma, 3)*(float(((int)(threadIdx.x / sz)) + 0.5f - y)*exp(-1.0f / 2.0f*pow((float(((int)(threadIdx.x / sz)) + 0.5f - y) / PSFsigma), 2.0f)) - float(((int)(threadIdx.x / sz)) - 0.5f - y)*exp(-1.0f / 2.0f*pow(float(((int)(threadIdx.x / sz)) - 0.5f - y) / PSFsigma, 2.0f)))*PSFx;
		deriv2[2] = 0.0f;
		deriv2[3] = 0.0f;

		//sync threads before sum
		__syncthreads();

		for (count1 = 0; count1 <vnum; count1++)
		{
			//clear s_arrx and s_arry to store for s_MM each pixel
			s_arrx[threadIdx.x] = 0.0;
			s_arry[threadIdx.x] = 0.0;
			__syncthreads();

			if ((model > 1e-3f) && (data > 0))
			{
				s_arrx[threadIdx.x] = deriv1[count1] * (data / model - 1);
				s_arry[threadIdx.x] = deriv2[count1] * (data / model - 1) - powf(deriv1[count1], 2)*(data / powf(model, 2));
			}
			__syncthreads();

			//sum up elements in s_arrx to calculate s_MM[count1] and sum up elements in s_arry to calculate s_MM[count1+5]
			for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
			{
				if (threadIdx.x < s)
				{
					s_arrx[threadIdx.x] += s_arrx[threadIdx.x + s];
					s_arry[threadIdx.x] += s_arry[threadIdx.x + s];
				}
				//sync threads before next loop
				__syncthreads();
				if (threadIdx.x == 0)
				{
					s_arrx[threadIdx.x] += s_arrx[2 * s];
					s_arry[threadIdx.x] += s_arry[2 * s];
					s_MM[count1] = s_arrx[threadIdx.x];
					s_MM[count1 + 5] = s_arry[threadIdx.x];
				}
			}
		}

		//sync threads before updating
		__syncthreads();

		rt = 1.0f;

		//update
		x -= rt*min(max(s_MM[0] / s_MM[0 + 5], -xlim), xlim);
		y -= rt*min(max(s_MM[1] / s_MM[1 + 5], -ylim), ylim);
		I -= rt*min(max(s_MM[2] / s_MM[2 + 5], -Ilim), Ilim);
		bg -= rt*min(max(s_MM[3] / s_MM[3 + 5], -bglim), bglim);
		I = max(I, 5.0f);
		bg = max(bg, 0.01f);

		//sync threads before next iteration
		__syncthreads();
	}


	//*******************************************************************************************
	//initialize
	if (threadIdx.x < 16) { s_MM[threadIdx.x] = 0.0; }
	s_arr[threadIdx.x] = 0.0;

	PSFx = 1.0f / 2.0f*(erf(float(((int)(threadIdx.x % sz)) - x + 0.5f)*sqrt(norm)) - erf(float(((int)(threadIdx.x % sz)) - x - 0.5f)*sqrt(norm)));
	PSFy = 1.0f / 2.0f*(erf(float(((int)(threadIdx.x / sz)) - y + 0.5f)*sqrt(norm)) - erf(float(((int)(threadIdx.x / sz)) - y - 0.5f)*sqrt(norm)));
	model = bg + I*PSFx * PSFy + r_varim / r_gainim / r_gainim;
	data = r_data + r_varim / r_gainim / r_gainim;

	//calculate derivatives
	deriv1[0] = -I / sqrt(2.0f*pi) / PSFsigma*(exp(-1.0f / 2.0f*pow((float(((int)(threadIdx.x % sz)) + 0.5f - x) / PSFsigma), 2.0f)) - expf(-1.0f / 2.0f*pow(float(((int)(threadIdx.x % sz)) - 0.5f - x) / PSFsigma, 2.0f)))*PSFy;
	deriv1[1] = -I / sqrt(2.0f*pi) / PSFsigma*(exp(-1.0f / 2.0f*pow((float(((int)(threadIdx.x / sz)) + 0.5f - y) / PSFsigma), 2.0f)) - exp(-1.0f / 2.0f*pow(float(((int)(threadIdx.x / sz)) - 0.5f - y) / PSFsigma, 2.0f)))*PSFx;
	deriv1[2] = PSFx * PSFy;
	deriv1[3] = 1.0f;

	//sync threads before building Fisher Information Matrix
	__syncthreads();

	//Building the Fisher Information Matrix
	for (countx = 0; countx < vnum; countx++)
	{
		for (county = countx; county < vnum; county++)
		{
			//sync threads to make sure threads are in same loop
			__syncthreads();
			s_arrx[threadIdx.x] = deriv1[county] * deriv1[countx] / model;
			//sync threads before sum
			__syncthreads();

			//calculate each element of s_MM. Use parallel reduction to speed up
			for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
			{
				if (threadIdx.x < s)
				{
					s_arrx[threadIdx.x] += s_arrx[threadIdx.x + s];
				}
				__syncthreads();
				if (threadIdx.x == 0)
				{
					s_arrx[threadIdx.x] += s_arrx[2 * s];
					s_MM[countx * 4 + county] = s_arrx[threadIdx.x];
					s_MM[county * 4 + countx] = s_MM[countx * 4 + county];
				}
			}
		}
	}

	if (model > 0)
	{
		if (data > 0){ s_arr[threadIdx.x] = data * logf(model) - model - data * logf(data) + data; }
		else{ s_arr[threadIdx.x] = -model; }
	}
	//sync threads before calculating Loglikelihood
	__syncthreads();

	//calculate Loglikelihood. Use parallel reduction to speed up
	for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (threadIdx.x < s)
		{
			s_arr[threadIdx.x] += s_arr[threadIdx.x + s];
		}
		//sync threads before next loop
		__syncthreads();
		if (threadIdx.x == 0)
		{
			s_arr[threadIdx.x] += s_arr[2 * s];
			llr = s_arr[threadIdx.x];
		}
	}

	//*******************************************************************************************

	if (threadIdx.x == 0)
	{
		//Matrix inverse
		kernel_MatInv(s_MM, s_invM, s_CRLB, vnum);
		//write to global memory
		d_x[blockIdx.x] = x;
		d_y[blockIdx.x] = y;
		d_I[blockIdx.x] = I;
		d_bg[blockIdx.x] = bg;
		d_ll[blockIdx.x] = llr;
	}

	//write to global memory
	if (threadIdx.x < vnum)
	{
		d_C[gridDim.x*threadIdx.x + blockIdx.x] = s_CRLB[threadIdx.x];
	}

}


__global__ void kernel_fit_two(float *d_data, float *d_varim, float *d_gainim,
	float PSFsigma, int sz, int iterations, float *d_x, float *d_y, float *d_I, float *d_bg, float *d_s, float *d_C, float *d_ll, int Nfits)
{
	//declare variables and set initial value
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int vnum = 5, count = 0, count1 = 0, countx = 0, county = 0, xx = 0, yy = 0;
	float rt = 0.0, PSFx = 0.0, PSFy = 0.0, model = 0.0, data = 0.0, r_data, r_varim, r_gainim, xlim = 1.0, ylim = 1.0, Ilim = 200.0, bglim = 2.0, slim = 0.3, deriv1[5] = { 0 }, deriv2[5] = { 0 }, s = 0.0, x, y, I = 0.0f, bg, tmp = 0.0f, tmpx = 0.0f, tmpy = 0.0f, sum = 0.0f, llr = 0.0f;

	__shared__ float s_MM[25], s_CRLB[5], s_invM[25];
	__shared__ float s_arr[441], s_arrx[441], s_arry[441];//Matrices to store temporal data before reduction

	//initialization of shared memory
	if (threadIdx.x < 25){ s_MM[threadIdx.x] = 0.0f; s_invM[threadIdx.x] = 0.0f; }
	if (threadIdx.x < 5){ s_CRLB[threadIdx.x] = 0.0f; }
	s_arr[threadIdx.x] = 0.0f;
	s_arrx[threadIdx.x] = 0.0f;
	s_arry[threadIdx.x] = 0.0f;

	//*******************************************************************************************
	//load data, variancemap and gainmap from global memory to register
	r_data = d_data[idx];
	r_varim = d_varim[idx];
	r_gainim = d_gainim[idx];

	//initial guess of x, y, background, photon and sigma
	s = PSFsigma;
	I = 0.0f;
	bg = 1e10;
	x = (float)sz / 2;
	y = (float)sz / 2;

	//loop through each pixel of one subregion to refine bg and photon initial guess. scmos noise ignored
	float norm = 1.0f / 2.0f / PSFsigma / PSFsigma;
	for (xx = 0; xx < sz; xx++)
	{
		for (yy = 0; yy < sz; yy++)
		{
			//initialize
			sum = 0;
			tmp = 0;

			//sync threads to make sure threads are in the same loop
			__syncthreads();

			//C++ use Row-major order. threadIdx%sz represents the column index here, which is row index in Matlab of a subregion
			s_arry[threadIdx.x] = exp(-pow((float)(((int)(threadIdx.x % sz)) - xx), 2)*norm)*exp(-pow((float)(yy - ((int)(threadIdx.x / sz))), 2)*norm);
			s_arrx[threadIdx.x] = s_arry[threadIdx.x] * r_data;

			//sync threads. So that s_arrx and s_arry is updated before sum
			__syncthreads();

			//calculate sum and tmp. Use parallel reduction algorithm to speed up
			for (unsigned int t = blockDim.x / 2; t > 0; t >>= 1)
			{
				if (threadIdx.x < t)
				{
					s_arrx[threadIdx.x] += s_arrx[threadIdx.x + t];
					s_arry[threadIdx.x] += s_arry[threadIdx.x + t];
				}
				//sync threads before next loop
				__syncthreads();

				if (threadIdx.x == 0)
				{
					s_arrx[threadIdx.x] += s_arrx[2 * t];
					s_arry[threadIdx.x] += s_arry[2 * t];
				}
			}

			//sync threads to wait for s_arrx[0] and s_arry[0] being calculated before load them to registers of each thread
			__syncthreads();

			tmp = s_arrx[0];
			sum = s_arry[0];
			tmp = tmp / sum;

			if (I < tmp)
			{
				x = (float)xx;
				y = (float)yy;
				I = tmp;
			}
			bg = min(bg, tmp);
		}
	}

	//constraint of bg and photon
	bg = max(bg, 2);
	I = max(0.0f, (I - bg) * 2 * pi*s*s);

	//*******************************************************************************************
	//iteration starts
	for (count = 0; count < iterations; count++)
	{
		//initialize
		norm = 1.0f / 2.0f / s / s;
		if (threadIdx.x < vnum + 5){ s_MM[threadIdx.x] = 0.0f; }
		s_arrx[threadIdx.x] = 0.0f;
		s_arry[threadIdx.x] = 0.0f;

		PSFx = 1.0f / 2.0f*(erf((((int)(threadIdx.x % sz)) - x + 0.5f)*sqrt(norm)) - erf((((int)(threadIdx.x % sz)) - x - 0.5f)*sqrt(norm)));
		PSFy = 1.0f / 2.0f*(erf(((int)(threadIdx.x / sz) - y + 0.5f)*sqrt(norm)) - erf(((int)(threadIdx.x / sz) - y - 0.5f)*sqrt(norm)));
		model = bg + I*PSFx*PSFy + r_varim / r_gainim / r_gainim;
		data = r_data + r_varim / r_gainim / r_gainim;

		//calculate derivatives
		deriv1[0] = -I / sqrt(2.0f*pi) / s*(exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x % sz)) + 0.5f - x) / s), 2.0f)) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x % sz)) - 0.5f - x) / s, 2.0f)))*PSFy;
		deriv1[1] = -I / sqrt(2.0f*pi) / s*(exp(-1.0f / 2.0f*pow((((int)(threadIdx.x / sz) + 0.5f - y) / s), 2.0f)) - exp(-1.0f / 2.0f*pow(((int)(threadIdx.x / sz) - 0.5f - y) / s, 2.0f)))*PSFx;
		deriv1[2] = PSFx*PSFy;
		deriv1[3] = 1.0f;
		tmpx = -I / sqrt(2.0f*pi) / s / s*(exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x % sz)) + 0.5f - x) / s), 2.0f))*(((int)(threadIdx.x % sz)) - x + 0.5f) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x % sz)) - 0.5f - x) / s, 2.0f))*(((int)(threadIdx.x % sz)) - x - 0.5f))*PSFy;
		tmpy = -I / sqrt(2.0f*pi) / s / s*(exp(-1.0f / 2.0f*pow((((int)(threadIdx.x / sz) + 0.5f - y) / s), 2.0f))*((int)(threadIdx.x / sz) - y + 0.5f) - exp(-1.0f / 2.0f*pow(((int)(threadIdx.x / sz) - 0.5f - y) / s, 2.0f))*((int)(threadIdx.x / sz) - y - 0.5f))*PSFx;
		deriv1[4] = 0;// tmpx + tmpy;

		deriv2[0] = -I / sqrt(2.0f*pi) / pow(s, 3)*((((int)(threadIdx.x % sz)) + 0.5f - x)*exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x % sz)) + 0.5f - x) / s), 2.0f)) - (((int)(threadIdx.x % sz)) - 0.5f - x)*exp(-1.0f / 2.0f*pow((((int)(threadIdx.x % sz)) - 0.5f - x) / s, 2.0f)))*PSFy;
		deriv2[1] = -I / sqrt(2.0f*pi) / pow(s, 3)*(((int)(threadIdx.x / sz) + 0.5f - y)*exp(-1.0f / 2.0f*pow((((int)(threadIdx.x / sz) + 0.5f - y) / s), 2.0f)) - ((int)(threadIdx.x / sz) - 0.5f - y)*exp(-1.0f / 2.0f*pow(((int)(threadIdx.x / sz) - 0.5f - y) / s, 2.0f)))*PSFx;
		deriv2[2] = 0.0f;
		deriv2[3] = 0.0f;
		deriv2[4] = -2.0f / s*tmpx - I / sqrt(2.0f*pi) / pow(s, 5)*(exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x % sz)) + 0.5f - x) / s), 2.0f))*pow((((int)(threadIdx.x % sz)) - x + 0.5f), 3) - exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x % sz)) - 0.5f - x) / s), 2.0f))*pow((((int)(threadIdx.x % sz)) - x - 0.5f), 3))*PSFy;
		deriv2[4] += -2.0f / s*tmpy - I / sqrt(2.0f*pi) / pow(s, 5)*(exp(-1.0f / 2.0f*pow((((int)(threadIdx.x / sz) + 0.5f - y) / s), 2.0f))*pow(((int)(threadIdx.x / sz) - y + 0.5f), 3) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x / sz) - 0.5f - y) / s), 2.0f))*pow(((int)(threadIdx.x / sz) - y - 0.5f), 3))*PSFx;

		//sync threads before sum
		__syncthreads();

		for (count1 = 0; count1 <vnum; count1++)
		{
			//initialize
			s_arrx[threadIdx.x] = 0.0;
			s_arry[threadIdx.x] = 0.0;
			__syncthreads();

			if ((model > 1e-3f) && (data > 0))
			{
				s_arrx[threadIdx.x] = deriv1[count1] * (data / model - 1);
				s_arry[threadIdx.x] = deriv2[count1] * (data / model - 1) - powf(deriv1[count1], 2)*(data / powf(model, 2));
			}
			__syncthreads();

			//calculate s_MM. Use parallel reduction to speed up
			for (unsigned int t = blockDim.x / 2; t > 0; t >>= 1)
			{
				if (threadIdx.x < t)
				{
					s_arrx[threadIdx.x] += s_arrx[threadIdx.x + t];
					s_arry[threadIdx.x] += s_arry[threadIdx.x + t];
				}
				//sync threads before next loop
				__syncthreads();

				if (threadIdx.x == 0)
				{
					s_arrx[threadIdx.x] += s_arrx[2 * t];
					s_arry[threadIdx.x] += s_arry[2 * t];
					s_MM[count1] = s_arrx[threadIdx.x];
					s_MM[count1 + 5] = s_arry[threadIdx.x];
				}
			}
		}

		//sync threads before update
		__syncthreads();

		rt = 1.0f;

		//update
		x -= rt*min(max(s_MM[0] / s_MM[0 + 5], -xlim), xlim);
		y -= rt*min(max(s_MM[1] / s_MM[1 + 5], -ylim), ylim);
		I -= rt*min(max(s_MM[2] / s_MM[2 + 5], -Ilim), Ilim);
		bg -= rt*min(max(s_MM[3] / s_MM[3 + 5], -bglim), bglim);
		s -= rt*min(max(s_MM[3] / s_MM[3 + 5], -slim), slim);

		I = max(I, 1.0f);
		bg = max(bg, 0.01f);
		s = min(s, sz / 1.5f);
		s = max(s, 0.3f);

		//sync threads before next loop
		__syncthreads();
	}

	//******************************************************************************************
	//initialize
	if (threadIdx.x < 25){ s_MM[threadIdx.x] = 0; }
	s_arr[threadIdx.x] = 0.0;

	PSFx = 1.0f / 2.0f*(erf((((int)(threadIdx.x % sz)) - x + 0.5f)*sqrt(norm)) - erf((((int)(threadIdx.x % sz)) - x - 0.5f)*sqrt(norm)));
	PSFy = 1.0f / 2.0f*(erf((((int)(threadIdx.x / sz)) - y + 0.5f)*sqrt(norm)) - erf((((int)(threadIdx.x / sz)) - y - 0.5f)*sqrt(norm)));
	model = bg + I*PSFx * PSFy + r_varim / r_gainim / r_gainim;
	data = r_data + r_varim / r_gainim / r_gainim;

	//calculate derivatives
	deriv1[0] = -I / sqrt(2.0f*pi) / s*(exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x % sz)) + 0.5f - x) / s), 2.0f)) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x % sz)) - 0.5f - x) / s, 2.0f)))*PSFy;
	deriv1[1] = -I / sqrt(2.0f*pi) / s*(exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x / sz)) + 0.5f - y) / s), 2.0f)) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x / sz)) - 0.5f - y) / s, 2.0f)))*PSFx;
	deriv1[2] = PSFx * PSFy;
	deriv1[3] = 1.0f;
	deriv1[4] = -I / sqrt(2.0f*pi) / s / s*(exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x % sz)) + 0.5f - x) / s), 2.0f))*(((int)(threadIdx.x % sz)) - x + 0.5f) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x % sz)) - 0.5f - x) / s, 2.0f))*(((int)(threadIdx.x % sz)) - x - 0.5f))*PSFy;
	deriv1[4] += -I / sqrt(2.0f*pi) / s / s*(exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x / sz)) + 0.5f - y) / s), 2.0f))*(((int)(threadIdx.x / sz)) - y + 0.5f) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x / sz)) - 0.5f - y) / s, 2.0f))*(((int)(threadIdx.x / sz)) - y - 0.5f))*PSFx;

	//snc threads before building the Fisher Information Matrix
	__syncthreads();

	//Building the Fisher Information Matrix
	for (countx = 0; countx < vnum; countx++)for (county = countx; county < vnum; county++)
	{
		//sync threads to make sure threads are in the same loop
		__syncthreads();
		s_arrx[threadIdx.x] = deriv1[county] * deriv1[countx] / model;
		//sync threads before sum
		__syncthreads();

		//calculate s_MM
		for (unsigned int t = blockDim.x / 2; t > 0; t >>= 1)
		{
			if (threadIdx.x < t)
			{
				s_arrx[threadIdx.x] += s_arrx[threadIdx.x + t];
			}
			__syncthreads();
			if (threadIdx.x == 0)
			{
				s_arrx[threadIdx.x] += s_arrx[2 * t];
				s_MM[countx * 5 + county] = s_arrx[threadIdx.x];
				s_MM[county * 5 + countx] = s_MM[countx * 5 + county];
			}
		}
	}

	if (model > 0)
	{
		if (data > 0) { s_arr[threadIdx.x] = data * log(model) - model - data * log(data) + data; }
		else{ s_arr[threadIdx.x] = -model; }
	}
	//snc threads before calculating Loglikelihood
	__syncthreads();

	//calculate Loglikelihood. Use parallel reduction to speed up
	for (unsigned int t = blockDim.x / 2; t > 0; t >>= 1)
	{
		if (threadIdx.x < t)
		{
			s_arr[threadIdx.x] += s_arr[threadIdx.x + t];
		}
		//sync threads before next loop
		__syncthreads();
		if (threadIdx.x == 0)
		{
			s_arr[threadIdx.x] += s_arr[2 * t];
			llr = s_arr[threadIdx.x];
		}
	}

	//*******************************************************************************************
	if (threadIdx.x == 0)
	{
		//Matrix inverse
		kernel_MatInv(s_MM, s_invM, s_CRLB, vnum);
		//write to global memory
		d_x[blockIdx.x] = x;
		d_y[blockIdx.x] = y;
		d_I[blockIdx.x] = I;
		d_bg[blockIdx.x] = bg;
		d_s[blockIdx.x] = s;
		d_ll[blockIdx.x] = llr;
	}
	//write to global memory
	if (threadIdx.x < vnum)
	{
		d_C[gridDim.x*threadIdx.x + blockIdx.x] = s_CRLB[threadIdx.x];
	}


}

__global__ void kernel_MLEFit_sigmaxy(float *d_data, float *d_varim, float *d_gainim, float PSFSigma, const int sz, const int iterations, float *d_C, float *d_ll, const int Nfits, float *d_x, float *d_y, float *d_I, float *d_bg, float *d_s, float *d_s2){
	//theta[6] is {x,y,photon,background,sigmax,sigmay}
	//declare variables and initial value
	float model = 0.0, cf = 0.0, df = 0.0, data = 0.0, Div = 0.0, PSFy = 0.0, PSFx = 0.0, dudt[6] = { 0.0 }, d2udt2[6] = { 0.0 }, maxjump[6] = { 1e0f, 1e0f, 1e2f, 2e0f, 1e-1f, 1e-1f }, g[6] = { 1.0f, 1.0f, 0.5f, 1.0f, 1.0f, 1.0f }, Nmax = 0.0, r_data = 0.0, r_varim = 0.0, r_gainim = 0.0, theta[6] = { 0.0 }, tmpx = 0.0f, tmpy = 0.0f, tmpsum = 0.0f;
	int NV = 6, idx = blockIdx.x*blockDim.x + threadIdx.x, count, count1, count2;
	int kk, ll;
	float filteredpixel = 0, sum = 0;

	__shared__ float s_arrx[441], s_arry[441], s_arr[441];
	__shared__ float s_NR_Numerator[6], s_NR_Denominator[6], s_M[6 * 6], s_Diag[6], s_Minv[6 * 6];

	//initialization of shared memory
	s_arrx[threadIdx.x] = 0;
	s_arry[threadIdx.x] = 0;
	s_arr[threadIdx.x] = 0;
	if (threadIdx.x < 36){ s_M[threadIdx.x] = 0; s_Minv[threadIdx.x] = 0; }

	//load data
	r_data = d_data[idx];
	r_varim = d_varim[idx];
	r_gainim = d_gainim[idx];

	//initial values
	s_arrx[threadIdx.x] = r_data * float((int)(threadIdx.x % sz));// elements for tmpx
	s_arry[threadIdx.x] = r_data * float((int)(threadIdx.x / sz));// elements for tmpy
	s_arr[threadIdx.x] = r_data;// elements for tmpsum

	//sync threads before sum
	__syncthreads();

	//calculate tmpx, tmpy, tmpsum.
	for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (threadIdx.x < s)
		{
			s_arrx[threadIdx.x] += s_arrx[threadIdx.x + s];
			s_arry[threadIdx.x] += s_arry[threadIdx.x + s];
			s_arr[threadIdx.x] += s_arr[threadIdx.x + s];
		}
		//sync threads before next loop
		__syncthreads();
		if (threadIdx.x == 0)
		{
			s_arrx[threadIdx.x] += s_arrx[2 * s];
			s_arry[threadIdx.x] += s_arry[2 * s];
			s_arr[threadIdx.x] += s_arr[2 * s];
		}
	}

	__syncthreads();

	tmpx = s_arrx[0];
	tmpy = s_arry[0];
	tmpsum = s_arr[0];

	theta[0] = tmpx / tmpsum;//x
	theta[1] = tmpy / tmpsum;//y

	//loop through all pixels to refine bg, photon, sigmax and sigmay initial guess. scmos noise ignored
	float norm = 1.0f / 2.0f / PSFSigma / PSFSigma;
	Nmax = 0.0f;
	theta[3] = 10e10f; //background
	for (kk = 0; kk < sz; kk++) for (ll = 0; ll < sz; ll++){

		//initialize
		filteredpixel = 0.0f;
		sum = 0.0f;

		//sync threads to make sure threads are in same loop
		__syncthreads();

		s_arry[threadIdx.x] = exp(-pow((float)(((int)(threadIdx.x % sz)) - kk), 2)*norm)*exp(-pow((float)(ll - ((int)(threadIdx.x / sz))), 2)*norm);//elements for sum
		s_arrx[threadIdx.x] = s_arry[threadIdx.x] * r_data;//elements for filteredpixel

		//sync threads before sum
		__syncthreads();

		for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
		{
			if (threadIdx.x < s)
			{
				s_arrx[threadIdx.x] += s_arrx[threadIdx.x + s];
				s_arry[threadIdx.x] += s_arry[threadIdx.x + s];
			}
			//sync threads before next loop
			__syncthreads();

			if (threadIdx.x == 0)
			{
				s_arrx[threadIdx.x] += s_arrx[2 * s];
				s_arry[threadIdx.x] += s_arry[2 * s];
			}
		}

		__syncthreads();

		filteredpixel = s_arrx[0];
		sum = s_arry[0];
		filteredpixel /= sum;

		Nmax = max(Nmax, filteredpixel);
		theta[3] = min(theta[3], filteredpixel);//background
	}
	theta[2] = max(0.0f, (Nmax - theta[3]) * 2 * pi*PSFSigma*PSFSigma);//photon
	theta[4] = PSFSigma;//sigmax
	theta[5] = PSFSigma;//sigmay

	//*******************************************************************************************************
	//iteration starts
	for (count = 0; count < iterations; count++)
	{
		//initialize
		if (threadIdx.x < 6){ s_NR_Numerator[threadIdx.x] = 0; s_NR_Denominator[threadIdx.x] = 0; }
		s_arrx[threadIdx.x] = 0.0;
		s_arry[threadIdx.x] = 0.0;

		norm = 1.0f / 2.0f / theta[4] / theta[4];
		PSFx = 1.0f / 2.0f*(erf((((int)(threadIdx.x % sz)) - theta[0] + 0.5f)*sqrt(norm)) - erf((((int)(threadIdx.x % sz)) - theta[0] - 0.5f)*sqrt(norm)));

		norm = 1.0f / 2.0f / theta[5] / theta[5];
		PSFy = 1.0f / 2.0f*(erf((((int)(threadIdx.x / sz)) - theta[1] + 0.5f)*sqrt(norm)) - erf((((int)(threadIdx.x / sz)) - theta[1] - 0.5f)*sqrt(norm)));

		model = theta[3] + theta[2] * PSFx*PSFy + r_varim / r_gainim / r_gainim;
		data = r_data + r_varim / r_gainim / r_gainim;

		//calculating derivatives
		dudt[0] = -theta[2] / sqrt(2.0f*pi) / theta[4] * (exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x % sz)) + 0.5f - theta[0]) / theta[4]), 2.0f)) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x % sz)) - 0.5f - theta[0]) / theta[4], 2.0f)))*PSFy;
		dudt[1] = -theta[2] / sqrt(2.0f*pi) / theta[5] * (exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x / sz)) + 0.5f - theta[1]) / theta[5]), 2.0f)) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x / sz)) - 0.5f - theta[1]) / theta[5], 2.0f)))*PSFx;

		dudt[2] = PSFx*PSFy;
		d2udt2[2] = 0.0f;
		dudt[3] = 1.0f;
		d2udt2[3] = 0.0f;

		dudt[4] = -theta[2] / sqrt(2.0f*pi) / theta[4] / theta[4] * (exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x % sz)) + 0.5f - theta[0]) / theta[4]), 2.0f))*(((int)(threadIdx.x % sz)) - theta[0] + 0.5f) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x % sz)) - 0.5f - theta[0]) / theta[4], 2.0f))*(((int)(threadIdx.x % sz)) - theta[0] - 0.5f))*PSFy;

		dudt[5] = -theta[2] / sqrt(2.0f*pi) / theta[5] / theta[5] * (exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x / sz)) + 0.5f - theta[1]) / theta[5]), 2.0f))*(((int)(threadIdx.x / sz)) - theta[1] + 0.5f) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x / sz)) - 0.5f - theta[1]) / theta[5], 2.0f))*(((int)(threadIdx.x / sz)) - theta[1] - 0.5f))*PSFx;

		if (model > 10e-3f) cf = data / model - 1;
		if (model > 10e-3f) df = data / pow(model, 2);
		cf = min(cf, 10e4f);
		df = min(df, 10e4f);

		for (count1 = 0; count1 < NV; count1++)
		{
			//sync threads to make sure threads are in same loop
			__syncthreads();

			s_arrx[threadIdx.x] = dudt[count1] * cf;//elements for NR_Numerator
			s_arry[threadIdx.x] = d2udt2[count1] * cf - pow(dudt[count1], 2)*df;//elements for NR_Denominator

			//syncthreads before sum
			__syncthreads();

			//calculate s_NR_Numerator and s_NR_Denominator
			for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
			{
				if (threadIdx.x < s)
				{
					s_arrx[threadIdx.x] += s_arrx[threadIdx.x + s];
					s_arry[threadIdx.x] += s_arry[threadIdx.x + s];
				}
				//sync threads before next loop 
				__syncthreads();
				if (threadIdx.x == 0)
				{
					s_arrx[threadIdx.x] += s_arrx[2 * s];
					s_arry[threadIdx.x] += s_arry[2 * s];
					s_NR_Numerator[count1] = s_arrx[threadIdx.x];
					s_NR_Denominator[count1] = s_arry[threadIdx.x];
				}
			}

		}

		//sync threads before update
		__syncthreads();

		// The update
		for (count2 = 0; count2 < NV; count2++){
			theta[count2] -= g[count2] * min(max(s_NR_Numerator[count2] / s_NR_Denominator[count2], -maxjump[count2]), maxjump[count2]);
		}

		// Any other constraints
		theta[2] = max(theta[2], 1.0f);//photon
		theta[3] = max(theta[3], 0.01f);//background
		theta[4] = max(theta[4], PSFSigma / 10.0f);//sigmax
		theta[5] = max(theta[5], PSFSigma / 10.0f);//sigmay

		//sync threads before next loop
		__syncthreads();
	}

	//*************************************************************************************************************
	// Calculating the CRLB and LogLikelihood
	//initialize
	Div = 0.0f;

	norm = 1.0f / 2.0f / theta[4] / theta[4];
	PSFx = 1.0f / 2.0f*(erf((((int)(threadIdx.x % sz)) - theta[0] + 0.5f)*sqrt(norm)) - erf((((int)(threadIdx.x % sz)) - theta[0] - 0.5f)*sqrt(norm)));

	norm = 1.0f / 2.0f / theta[5] / theta[5];
	PSFy = 1.0f / 2.0f*(erf((((int)(threadIdx.x / sz)) - theta[1] + 0.5f)*sqrt(norm)) - erf((((int)(threadIdx.x / sz)) - theta[1] - 0.5f)*sqrt(norm)));
	model = theta[3] + theta[2] * PSFx*PSFy + r_varim / r_gainim / r_gainim;
	data = r_data + r_varim / r_gainim / r_gainim;

	//calculating derivatives
	dudt[0] = -theta[2] / sqrt(2.0f*pi) / theta[4] * (exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x % sz)) + 0.5f - theta[0]) / theta[4]), 2.0f)) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x % sz)) - 0.5f - theta[0]) / theta[4], 2.0f)))*PSFy;
	dudt[1] = -theta[2] / sqrt(2.0f*pi) / theta[5] * (exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x / sz)) + 0.5f - theta[1]) / theta[5]), 2.0f)) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x / sz)) - 0.5f - theta[1]) / theta[5], 2.0f)))*PSFx;
	dudt[2] = PSFx*PSFy;
	dudt[3] = 1.0f;
	dudt[4] = -theta[2] / sqrt(2.0f*pi) / theta[4] / theta[4] * (exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x % sz)) + 0.5f - theta[0]) / theta[4]), 2.0f))*(((int)(threadIdx.x % sz)) - theta[0] + 0.5f) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x % sz)) - 0.5f - theta[0]) / theta[4], 2.0f))*(((int)(threadIdx.x % sz)) - theta[0] - 0.5f))*PSFy;
	dudt[5] = -theta[2] / sqrt(2.0f*pi) / theta[5] / theta[5] * (exp(-1.0f / 2.0f*pow(((((int)(threadIdx.x / sz)) + 0.5f - theta[1]) / theta[5]), 2.0f))*(((int)(threadIdx.x / sz)) - theta[1] + 0.5f) - exp(-1.0f / 2.0f*pow((((int)(threadIdx.x / sz)) - 0.5f - theta[1]) / theta[5], 2.0f))*(((int)(threadIdx.x / sz)) - theta[1] - 0.5f))*PSFx;

	//sync threads before building the Fisher Information Matrix
	__syncthreads();

	//Building the Fisher Information Matrix
	for (count1 = 0; count1 < NV; count1++)for (count2 = count1; count2 < NV; count2++)
	{
		//sync threads to make sure threads are in same loop
		__syncthreads();
		s_arrx[threadIdx.x] = 0;
		s_arrx[threadIdx.x] = dudt[count2] * dudt[count1] / model;
		//sync threads before sum
		__syncthreads();

		//calculate s_M
		for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
		{
			if (threadIdx.x < s)
			{
				s_arrx[threadIdx.x] += s_arrx[threadIdx.x + s];
			}
			__syncthreads();
			if (threadIdx.x == 0)
			{
				s_arrx[threadIdx.x] += s_arrx[2 * s];
				s_M[count1 * NV + count2] = s_arrx[threadIdx.x];
				s_M[count2 * NV + count1] = s_M[count1 * NV + count2];
			}
		}
	}

	s_arr[threadIdx.x] = 0.0;
	if (model > 0)
	{
		if (data > 0){ s_arr[threadIdx.x] = data*log(model) - model - data*log(data) + data; }
		else{ s_arr[threadIdx.x] = -model; }
	}
	//sync threads before sum
	__syncthreads();

	//calculate LogLikelyhood
	for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
	{
		if (threadIdx.x < s)
		{
			s_arr[threadIdx.x] += s_arr[threadIdx.x + s];
		}
		//sync threads before next loop
		__syncthreads();
		if (threadIdx.x == 0)
		{
			s_arr[threadIdx.x] += s_arr[2 * s];
			Div = s_arr[threadIdx.x];
		}
	}

	if (threadIdx.x == 0)
	{
		// Matrix inverse
		kernel_MatInv(s_M, s_Minv, s_Diag, NV);
		//write to global memory
		d_ll[blockIdx.x] = Div;
		d_x[blockIdx.x] = theta[0];
		d_y[blockIdx.x] = theta[1];
		d_I[blockIdx.x] = theta[2];
		d_bg[blockIdx.x] = theta[3];
		d_s[blockIdx.x] = theta[4];
		d_s2[blockIdx.x] = theta[5];

	}
	//write to global memory
	if (threadIdx.x < NV)
	{
		d_C[gridDim.x*threadIdx.x + blockIdx.x] = s_Diag[threadIdx.x];
	}
	return;
}

