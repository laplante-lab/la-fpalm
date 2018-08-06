#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "cHistRecon.h"

//#include "matrix.h"

// Thread block size
#define CBLK 4000
#define pi 3.141592f

#define max(a,b) ( (a) >= (b) ? (a) : (b) )  
#define min(a,b) ( (a) < (b) ? (a) : (b) )  

uint16_t* cHistRecon(int sizeX, int sizeY, int size, float* X, float* Y, int useHE) {
  uint16_t *out;
  float *out2;
  int i;
  int N;
  //const mwSize *datasize;  
  int Xtmp,Ytmp;
  //mwSize outsize[2];
  int *Hist;
  float *CDF;
  int sum;
  //int cdf_min;
  int Npixels;
  const int NBIN = 1023;

  //get variables
  //datasize=mxGetDimensions(prhs[2]);

  //create output

  //outsize[0]=sizeX;
  //outsize[1]=sizeY;
  //plhs[0]=mxCreateNumericArray(2,outsize,mxUINT16_CLASS,mxREAL);
  //plhs[1]=mxCreateNumericArray(2,outsize,mxSINGLE_CLASS,mxREAL);


  out = (uint16_t *) calloc ( 2*sizeX*sizeY, sizeof ( uint16_t ) ); 
  out2 = (float *) calloc ( 2*sizeX*sizeY, sizeof(float) );
  // out = (uint16_t *) malloc ( sizeX*sizeY * sizeof ( uint16_t ) ); 
  //out2 = (float *) malloc ( sizeX*sizeY * sizeof(float) );
  //out=(uint16_T *) mxGetData(plhs[0]);
  //out2= (float *) mxGetData(plhs[1]);

  sum=0;
  //printf("size: %d\n", size);
  for (i=0;i<size;i++){
    Xtmp=max(0,min(sizeX-1,(int)floor(X[i])));
    Ytmp=max(0,min(sizeY-1,(int)floor(Y[i])));
//    printf("i(%d): Xtmp:%d Ytmp:%d\n", i, Xtmp, Ytmp);

    out[Ytmp*sizeX+Xtmp]+=1;
  }

  //printf("check1\n");

  if (useHE){ //Histogram equalization
    //Hist = (int*) mxCalloc(NBIN,sizeof(int));
    //CDF = (float*) mxCalloc(NBIN,sizeof(float));

    //printf("check2\n");
    Hist = (int*) calloc( NBIN, sizeof(int) );
    CDF = (float*) calloc( NBIN, sizeof(float) );
    for (i=0; i<NBIN; i++) Hist[i]=0;

    //printf("check3\n");
    Npixels=sizeX*sizeY;
    for (i=0;i<Npixels;i++) if(out[i]) Hist[out[i]]+=1; //calc histogram

    i=1; //ignoring zero pixels
    CDF[0] = 0;
    sum=0;
    //printf("check4\n");
    while (sum<size) { //calc CDF
      CDF[i]=CDF[i-1]+(float)Hist[i];
      sum+=Hist[i]*i;
      i++;
    }
    N = (int) CDF[i-1];
    //printf("%f %d %d\n",CDF[i-1],i-1,N);
    //printf("%f %f %f %f %f %f\n",CDF[0],CDF[1],CDF[2],CDF[3],CDF[4],CDF[5]);

    //out = (uint16_t*) realloc( out, sizeX*sizeY*sizeof(uint16_t) );
    //out2 = (float*) realloc( out2, sizeX*sizeY*sizeof(float) );
    //out = (uint16_t*) calloc(sizeX*sizeY, sizeof (uint16_t)); 
    //out2 = (float*) calloc(sizeX*sizeY, sizeof(float));

    //convert pixels
    for (i=0;i<Npixels;i++) if(out[i]) {
      out2[i]= NBIN*(CDF[out[i]]-1)/(N-1);
      out[i]=  (uint16_t) NBIN*(CDF[out[i]]-1)/(N-1);

      //printf("%0.4lf \t %d\n", out2[i], out[i]);
    }

    free(Hist);
    free(CDF);
  }

  free(out2);
  return out;
}
