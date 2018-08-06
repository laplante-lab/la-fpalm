/*
 * \File : definitions.h
 * \Author : Peiyi Zhang
 * \Date : November 3rd, 2016
 * \This is a modified version of Fang Huang's definitions.h in sCMOS software. 
 * \Copyright @ Purdue University
 */
#ifndef DEFINITIONS_H
#define DEFINITIONS_H
#include <cuda_runtime.h>

#define pi 3.141592f
void CUDAERROR(const char *instr,int lineNumber);
void cudasafe( cudaError_t err, char* str, int lineNumber);
#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#endif