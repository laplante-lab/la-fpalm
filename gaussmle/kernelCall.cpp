#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cuda_runtime.h>
#include "definitions.h"

#include "kernelCall.h"
#include "kernel.h"

//*******************************************************************************************
void CUDAERROR(const char *instr, int lineNumber)
{
  cudaError_t errornum;
  const char *str;
  if (errornum = cudaGetLastError()) {
    //reset all cuda devices
    int deviceCount = 0;
    int ii = 0;
    cudasafe(cudaGetDeviceCount(&deviceCount), (char *)"cudaGetDeviceCount", __LINE__); //query number of GPUs
    for (ii = 0; ii < deviceCount; ii++) {
      cudaSetDevice(ii);
      //cudaDeviceReset();
    }
    str = cudaGetErrorString(errornum);
    //cudaDeviceReset();
    //mexErrMsgIdAndTxt("mexFunction:cudaFail", "SRsCMOS_MLE(line %i): %s in %s\n", lineNumber, str, instr);
    exit(1); // might not stop matlab
  }
}

void cudasafe(cudaError_t err, char* str, int lineNumber)
{
  //printf("%s\n", str);
  if (err != cudaSuccess)
  {
    printf("FAILED\n");
    //reset all cuda devices
    int deviceCount = 0;
    int ii = 0;
    cudasafe(cudaGetDeviceCount(&deviceCount), (char *)"cudaGetDeviceCount", __LINE__); //query number of GPUs
    //for (ii = 0; ii < deviceCount; ii++) {
    //  cudaSetDevice(ii);
    //  cudaDeviceReset();
    //}
    printf("cudaFail %s failed with error code %i at line %d\n", str, err, lineNumber);
    //mexErrMsgIdAndTxt("mexFunction:cudaFail", "%s failed with error code %i at line %d\n", str, err, lineNumber);
    exit(1); // might not stop matlab
  }
}

extern "C" void kernel_fit_one_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float *d_varim, float *d_gainim, float gsig, int subsz, int iter_num,
float *d_x, float *d_y, float *d_I, float *d_bg, float *d_C, float *d_ll, int Nfits);

extern "C" void kernel_fit_two_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float *d_varim, float *d_gainim, float gsig, int subsz, int iter_num,
float *d_x, float *d_y, float *d_I, float *d_bg, float *d_s, float *d_C, float *d_ll, int Nfits);

extern void kernel_MLEFit_sigmaxy_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float *d_varim, float *d_gainim, float PSFSigma, const int sz, const int iterations, float *d_C, float *d_ll, const int Nfits, float *d_x, float *d_y, float *d_I, float *d_bg, float *d_s, float *d_s2);


void SRsCMOS_MLE( int num_dataset, int subtot, float *subregion, float *subvarim, float *subgainim, int fit_type, float gsig, int iter_num, float *sub_x, float *sub_y, float *sub_photon, float *sub_background, float *CRLB2, float *LL2 ) {

  //query GPUs
  int deviceCount = 0, driverVersion = 0, runtimeVersion = 0;
  cudasafe(cudaDriverGetVersion(&driverVersion), (char *)"Could not query CUDA driver version", __LINE__);
  cudasafe(cudaRuntimeGetVersion(&runtimeVersion), (char *)"Could not query CUDA runtime version", __LINE__);
  cudasafe(cudaGetDeviceCount(&deviceCount), (char *)"Error detecting CUDA devices", __LINE__);

  if (deviceCount < 1)
    printf("No CUDA capable devices were detected\n");

  //printf("Number of GPUs: %d\n", deviceCount);

  cudaDeviceProp deviceProp;

  //cudasafe(cudaSetDevice(0), (char *)"Could not select GPU 0.", __LINE__);
  //printf("Using GPU %s\n", deviceProp.name);
  //if (deviceProp.kernelExecTimeoutEnabled)
  //  printf("Warning, Kernel Execution Timeout is enabled for the GPU you are using.\nIf your fitting takes longer than the timeout it will fail.\n");

  if ( (deviceCount > 1) && (num_dataset % 2 == 0) ) {
    //printf("using second GPU\n");
    cudasafe(cudaSetDevice(1), (char *)"Could not select GPU 1.", __LINE__);
    cudasafe(cudaGetDeviceProperties(&deviceProp, 1), (char *)"Could not get properties for device 1.", __LINE__);
  }
  else {
    cudasafe(cudaSetDevice(0), (char *)"Could not select GPU 0.", __LINE__);
    cudasafe(cudaGetDeviceProperties(&deviceProp, 0), (char *)"Could not get properties for device 0.", __LINE__);
  }

  //declare all vars
  int Nfits, flag, gridsz, blocksz;
  //const mwSize *datasz;
  //float gsig;
  float *data, *d_data, *varim, *d_varim, *gainim, *d_gainim, *d_x, *d_y, *d_I, *d_bg, *d_s, *d_s2, *d_C, *d_ll;
  size_t subsz;

  //size_t ndim, subtot, subsz;
  //Using stream may help speed up the program a little bit. Because different streams may execute their commands out of order with respect to one another or concurrently. But the behavior is not guaranteed.
  //printf("creating streams\n");

  cudaStream_t stream0, stream1, stream2, stream3, stream4, stream5, stream6, stream7;
  cudaStreamCreate(&stream0);
  cudaStreamCreate(&stream1);
  cudaStreamCreate(&stream2);
  cudaStreamCreate(&stream3);
  cudaStreamCreate(&stream4);
  cudaStreamCreate(&stream5);
  cudaStreamCreate(&stream6);
  cudaStreamCreate(&stream7);

  //printf("DONE STREAMS\n");

  //datasz = mxGetDimensions(prhs[0]);
  //ndim = 2;
  //ndim = mxGetNumberOfDimensions(prhs[0]);

  //subtot = datasz[2];//number of subregions
  subsz = 7;//the sizes of the first dimension
  //subsz = datasz[0];//the sizes of the first dimension

  int BlockSize = subsz*subsz;

  Nfits = (int)subtot;

  data = subregion;
  varim = subvarim;
  gainim = subgainim;
  flag = fit_type;
  gsig = gsig;
  iter_num = iter_num;

  // data = (float *)mxGetData(prhs[0]);
  // varim = (float *)mxGetData(prhs[1]);
  // gainim = (float *)mxGetData(prhs[2]);
  // flag = (int)mxGetScalar(prhs[3]);
  // gsig = (float)mxGetScalar(prhs[4]);
  // iter_num = (int)mxGetScalar(prhs[5]);

  //allocate memory on the GPU, set initial value and copy variables from CPU to GPU
  //printf("do: %d\n", num_dataset);
  cudasafe(cudaMalloc((void**)&d_data, subsz*subsz*Nfits*sizeof(float)), (char *)"malloc d_data", __LINE__);
  //printf("done: %d\n", num_dataset);
  cudasafe(cudaMemsetAsync(d_data, 0, subsz*subsz*Nfits*sizeof(float), stream0), (char *)"memset d_data", __LINE__);
  cudasafe(cudaMemcpyAsync(d_data, data, subsz*subsz*Nfits*sizeof(float), cudaMemcpyHostToDevice, stream0), (char *)"memcopy d_data from data", __LINE__);

  cudasafe(cudaMalloc((void**)&d_varim, subsz*subsz*Nfits*sizeof(float)), (char *)"malloc d_varim", __LINE__);
  cudasafe(cudaMemsetAsync(d_varim, 0, subsz*subsz*Nfits*sizeof(float), stream1), (char *)"memset d_varim", __LINE__);
  cudasafe(cudaMemcpyAsync(d_varim, varim, subsz*subsz*Nfits*sizeof(float), cudaMemcpyHostToDevice, stream1), (char *)"memcopy d_varim from varim", __LINE__);

  cudasafe(cudaMalloc((void**)&d_gainim, subsz*subsz*Nfits*sizeof(float)), (char *)"malloc d_gainim", __LINE__);
  cudasafe(cudaMemsetAsync(d_gainim, 0, subsz*subsz*Nfits*sizeof(float), stream2), (char *)"memset d_gainim", __LINE__);
  cudasafe(cudaMemcpyAsync(d_gainim, gainim, subsz*subsz*Nfits*sizeof(float), cudaMemcpyHostToDevice, stream2), (char *)"memcopy d_gainim from gainim", __LINE__);

  // configure kernel
  gridsz = (int)subtot;
  blocksz = BlockSize;
  dim3 dimBlock(blocksz);
  dim3 dimGrid(gridsz);

  //if ((flag != 1) && (flag != 2) && (flag != 3))
  //  mexErrMsgTxt("Please choose between fitting type 1 , 2 and 3.");

  size_t cudabytes = (size_t)Nfits*sizeof(float);
  int outsz = Nfits;//subregion numbers

//mwSize outsz = (mwSize)Nfits;//subregion numbers
  if (flag == 1) {
    //printf("num of threads per block: %d,num of blocks: %d, Nfits: %d, subtot: %d\n", blocksz, gridsz, Nfits, subtot);

    //allocate memory on the GPU
    cudasafe(cudaMalloc((void**)&d_x, cudabytes), (char *)"malloc d_x", __LINE__);
    cudasafe(cudaMalloc((void**)&d_y, cudabytes), (char *)"malloc d_y", __LINE__);
    cudasafe(cudaMalloc((void**)&d_I, cudabytes), (char *)"malloc d_I", __LINE__);
    cudasafe(cudaMalloc((void**)&d_bg, cudabytes), (char *)"malloc d_bg", __LINE__);
    cudasafe(cudaMalloc((void**)&d_C, 4 * cudabytes), (char *)"malloc d_C", __LINE__);
    cudasafe(cudaMalloc((void**)&d_ll, cudabytes), (char *)"malloc d_ll", __LINE__);

    //create output matrix
    // plhs[0] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
    // plhs[1] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
    // plhs[2] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
    // plhs[3] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
    // plhs[4] = mxCreateNumericMatrix(outsz, 4, mxSINGLE_CLASS, mxREAL);
    // plhs[5] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
    //sub_x = (float *)malloc( subtot * sizeof (float) );
    //sub_y = (float *)malloc( subtot * sizeof (float) );
    //sub_photon = (float *)malloc( subtot * sizeof (float) );
    //sub_background = (float *)malloc( subtot * sizeof (float) );
    //CRLB2 = (float *)malloc( subtot * sizeof (float) );
    //LL2 = (float *)malloc( subtot * sizeof (float) );
 
    //launch kernel
    kernel_fit_one_wrapper(dimGrid, dimBlock, d_data, d_varim, d_gainim, gsig, (int)subsz, iter_num, d_x, d_y, d_I, d_bg, d_C, d_ll, Nfits);
    
    CUDAERROR("kernel_fit_one", __LINE__);

    //copy back from GPU to CPU
    cudasafe(cudaMemcpyAsync(sub_x, d_x, cudabytes, cudaMemcpyDeviceToHost, stream0), (char *)"copy d_x to matlab", __LINE__);
    cudasafe(cudaMemcpyAsync(sub_y, d_y, cudabytes, cudaMemcpyDeviceToHost, stream1), (char *)"copy d_y to matlab", __LINE__);
    cudasafe(cudaMemcpyAsync(sub_photon, d_I, cudabytes, cudaMemcpyDeviceToHost, stream2), (char *)"copy d_I to matlab", __LINE__);
    cudasafe(cudaMemcpyAsync(sub_background, d_bg, cudabytes, cudaMemcpyDeviceToHost, stream3), (char *)"copy d_bg to matlab", __LINE__);
    cudasafe(cudaMemcpyAsync(CRLB2, d_C, 4 * cudabytes, cudaMemcpyDeviceToHost, stream4), (char *)"copy d_C to matlab", __LINE__);
    cudasafe(cudaMemcpyAsync(LL2, d_ll, cudabytes, cudaMemcpyDeviceToHost, stream5), (char *)"copy d_ll to matlab", __LINE__);
  }

//if (flag == 2)
//{
//printf("num of threads per block: %d,num of blocks: %d, Nfits: %d, subtot: %d\n", blocksz, gridsz, Nfits, subtot);
//
////allocate memory on the GPU
//cudasafe(cudaMalloc((void**)&d_x, cudabytes), "malloc d_x", __LINE__);
//cudasafe(cudaMalloc((void**)&d_y, cudabytes), "malloc d_y", __LINE__);
//cudasafe(cudaMalloc((void**)&d_I, cudabytes), "malloc d_I", __LINE__);
//cudasafe(cudaMalloc((void**)&d_bg, cudabytes), "malloc d_bg", __LINE__);
//cudasafe(cudaMalloc((void**)&d_s, cudabytes), "malloc d_s", __LINE__);
//cudasafe(cudaMalloc((void**)&d_C, 5 * cudabytes), "malloc d_C", __LINE__);
//cudasafe(cudaMalloc((void**)&d_ll, cudabytes), "malloc d_ll", __LINE__);
//
////create output matrix
//// plhs[0] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//// plhs[1] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//// plhs[2] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//// plhs[3] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//// plhs[4] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//// plhs[5] = mxCreateNumericMatrix(outsz, 5, mxSINGLE_CLASS, mxREAL);
////plhs[6] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//
////launch kernel
//kernel_fit_two_wrapper(dimGrid, dimBlock, d_data, d_varim, d_gainim, gsig, (int)subsz, iter_num, d_x, d_y, d_I, d_bg, d_s, d_C, d_ll, Nfits);
//CUDAERROR("kernel_fit_two", __LINE__);
//
////copy back from GPU to CPU
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[0]), d_x, cudabytes, cudaMemcpyDeviceToHost, stream0), "copy d_x to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[1]), d_y, cudabytes, cudaMemcpyDeviceToHost, stream1), "copy d_y to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[2]), d_I, cudabytes, cudaMemcpyDeviceToHost, stream2), "copy d_I to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[3]), d_bg, cudabytes, cudaMemcpyDeviceToHost, stream3), "copy d_bg to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[4]), d_s, cudabytes, cudaMemcpyDeviceToHost, stream4), "copy d_s to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[5]), d_C, 5 * cudabytes, cudaMemcpyDeviceToHost, stream5), "copy d_C to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[6]), d_ll, cudabytes, cudaMemcpyDeviceToHost, stream6), "copy d_ll to matlab", __LINE__);
//
////clear up!
//cudasafe(cudaFree(d_s), "freeing d_s", __LINE__);
//}
//
//if (flag == 3)
//{
//printf("num of threads per block: %d,num of blocks: %d, Nfits: %d, subtot: %d\n", blocksz, gridsz, Nfits, subtot);
//
////allocate memory on the GPU
//cudasafe(cudaMalloc((void**)&d_C, 6 * cudabytes), "Failed cudaMalloc on d_C.", __LINE__);
//cudasafe(cudaMemset(d_C, 0, 6 * cudabytes), "Failed cudaMemset on d_C.", __LINE__);
//cudasafe(cudaMalloc((void**)&d_ll, cudabytes), "Failed cudaMalloc on d_ll.", __LINE__);
//cudasafe(cudaMalloc((void**)&d_x, cudabytes), "malloc d_x", __LINE__);
//cudasafe(cudaMalloc((void**)&d_y, cudabytes), "malloc d_y", __LINE__);
//cudasafe(cudaMalloc((void**)&d_I, cudabytes), "malloc d_I", __LINE__);
//cudasafe(cudaMalloc((void**)&d_bg, cudabytes), "malloc d_bg", __LINE__);
//cudasafe(cudaMalloc((void**)&d_s, cudabytes), "malloc d_s", __LINE__);
//cudasafe(cudaMalloc((void**)&d_s2, cudabytes), "malloc d_s2", __LINE__);
//
////create output matrix
//// plhs[0] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//// plhs[1] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//// plhs[2] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//// plhs[3] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//// plhs[4] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//// plhs[5] = mxCreateNumericMatrix(outsz, 1, mxSINGLE_CLASS, mxREAL);
//// plhs[6] = mxCreateNumericMatrix(Nfits, 6, mxSINGLE_CLASS, mxREAL);
//// plhs[7] = mxCreateNumericMatrix(Nfits, 1, mxSINGLE_CLASS, mxREAL);
//
////launch kernel
//kernel_MLEFit_sigmaxy_wrapper(dimGrid, dimBlock, d_data, d_varim, d_gainim, gsig, subsz, iter_num, d_C, d_ll, Nfits, d_x, d_y, d_I, d_bg, d_s, d_s2);
//CUDAERROR("kernel_fit_three", __LINE__);
//
////copy back from GPU to CPU
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[0]), d_x, cudabytes, cudaMemcpyDeviceToHost, stream0), "copy d_x to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[1]), d_y, cudabytes, cudaMemcpyDeviceToHost, stream1), "copy d_y to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[2]), d_I, cudabytes, cudaMemcpyDeviceToHost, stream2), "copy d_I to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[3]), d_bg, cudabytes, cudaMemcpyDeviceToHost, stream3), "copy d_bg to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[4]), d_s, cudabytes, cudaMemcpyDeviceToHost, stream4), "copy d_s to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((float *)mxGetData(plhs[5]), d_s2, cudabytes, cudaMemcpyDeviceToHost, stream5), "copy d_s2 to matlab", __LINE__);
//cudasafe(cudaMemcpyAsync((void*)mxGetData(plhs[6]), d_C, 6 * cudabytes, cudaMemcpyDeviceToHost, stream6), "cudaMemcpy failed for d_C.", __LINE__);
//cudasafe(cudaMemcpyAsync((void*)mxGetData(plhs[7]), d_ll, cudabytes, cudaMemcpyDeviceToHost, stream7), "cudaMemcpy failed for d_ll.", __LINE__);
//
////clean up!
//cudasafe(cudaFree(d_s), "freeing d_s", __LINE__);
//cudasafe(cudaFree(d_s2), "freeing d_s2", __LINE__);
//}

  //clean up!
  //printf("cuda clean up\n");
  cudasafe(cudaFree(d_varim), (char *)"freeing d_varim", __LINE__);
  cudasafe(cudaFree(d_gainim), (char *)"freeing d_gainim", __LINE__);
  cudasafe(cudaFree(d_data), (char *)"freeing d_data", __LINE__);
  cudasafe(cudaFree(d_x), (char *)"freeing d_x", __LINE__);
  cudasafe(cudaFree(d_y), (char *)"freeing d_y", __LINE__);
  cudasafe(cudaFree(d_I), (char *)"freeing d_I", __LINE__);
  cudasafe(cudaFree(d_bg), (char *)"freeing d_bg", __LINE__);
  cudasafe(cudaFree(d_ll), (char *)"freeing d_ll", __LINE__);
  cudasafe(cudaFree(d_C), (char *)"freeing d_C", __LINE__);

  cudaStreamDestroy(stream0);
  cudaStreamDestroy(stream1);
  cudaStreamDestroy(stream2);
  cudaStreamDestroy(stream3);
  cudaStreamDestroy(stream4);
  cudaStreamDestroy(stream5);
  cudaStreamDestroy(stream6);
  cudaStreamDestroy(stream7);

  //cudaDeviceReset();

  return;
}
