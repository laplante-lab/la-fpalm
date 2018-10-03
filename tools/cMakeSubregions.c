#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "cMakeSubregions.h"

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

// Thread block size

void getbox(const float *data,const int ii,const int sz,const int szinX,const int szinY,
        const int szinZ,const double *X,const double *Y,const double *Z, float* dataout,float* leftcoord,float* topcoord);
void version();

//void mexFunction(int nlhs, mxArray *plhs[],	int	nrhs, const	mxArray	*prhs[])
/** \brief This function takes a list of x,y,z coordinates and crops out subregions of size 
  * sz x sz at the specified locations in the dataset.  It returns the array of subregions along with a list
  * of left and top coordinates for each. 
  *  \param nlhs number of left hand mxArrays to return
  *  \param plhs array of pointers to the output mxArrays
  *  \param nrhs number of input mxArrays
  *  \param prhs array of pointers to the input mxArrays.
  */
float* cMakeSubregions( double *x, double *y, double *z, int num_regions, int sz, float* input_data, float* l2, float* t2, int szinX, int szinY, int szinZ ) {

  float *dataout = (float*)malloc(sz*sz*num_regions*sizeof(float));				//!< an array of 2D subregions
  float *leftcoord = t2; //!< an array of the left coordinates of the subregions
  float *topcoord = l2; //!< an array of the top coordinates of the subregions
 
  int ii;						//!< loop counter
  double *X = x;			//!< array of x coordinates
  double *Y = y;			//!< array of y coordinates
  double *Z = z;			//!< array of z coordinates

  float *data = input_data;			//!< dataset to crop regions out of

  int N = num_regions;

  for (ii=0;ii<N;ii++)
    getbox(data,ii,sz,szinX,szinY,szinZ,X,Y,Z,dataout,leftcoord,topcoord);

  return dataout; 
};

/** \brief This function copies the specified subregion in data to data 
  *  \param data the original data set to crop from
  *  \param ii the point to copy
  *  \param sz the size of the subregion to copy
  *  \param szinX x dimension of data
  *  \param szinY y dimension of data
  *  \param szinZ z dimension of data
  *  \param X x coordinate of the center of the subregion
  *  \param Y y coordinate of the center of the subregion
  *  \param Z z coordinate of the center of the subregion
  *  \param dataout array of subregions to copy to
  *  \param leftcoord left coordinate of the subregion in the original image
  *  \param topcoord right coordinate of the subregion in the original image
  */
void getbox(const float *data,const int ii,const int sz,const int szinX,const int szinY,
const int szinZ,const double *X,const double *Y,const double *Z, float* dataout,float* leftcoord,float* topcoord) {
//  int xx,yy;
  int yy;
  int l,r,t,b;

  const int szl=(int) floor(sz/2.0+.5);
  //get coordinates
  const int x=(int) floor(X[ii]+.5)+1; 
  const int y=(int) floor(Y[ii]+.5)+1;
  const int z=(int) Z[ii];

//#ifdef _DEBUG
  // x and y get +1 so they can = szin
  if (x<1 || y<1 || z<0 || x>szinX || y>szinY || z>szinZ) {
    printf("Point %d out of bounds position %d,%d,%d dataset size %d,%d,%d\n",ii,x,y,z,szinX, szinY, szinZ);
    //static char msgbuf[256];
    //sprintf(msgbuf,"Point %d out of bounds position %d,%d,%d dataset size %d,%d,%d",ii,x,y,z,szinX, szinY, szinZ);
    //mexErrMsgIdAndTxt("cMakeSubregions:InputOutofRange",msgbuf);
    exit(EXIT_FAILURE);
  }
//#endif

  //mexPrintf("ii: %d %d %d %d %d\n",ii,x,y,z,szl);
  //printf("ii: %d %d %d %d %d\n",ii,x,y,z,szl);

  //if (z>szinZ-1)
  //{
  //	mexPrintf("%d %d %d %d\n",x,y,z,szl);
  //    mexErrMsgTxt(" Z index larger than datasize\n");
  //}

  l=max(x-szl,0);
  r=l+sz-1;
  if (r>(szinX-1)) {r=szinX-1;l=r-sz+1;}
  t=max(y-szl,0);
  b=t+sz-1;
  if (b>(szinY-1)) {b=szinY-1;t=b-sz+1;}

  for (yy=0;yy<sz;yy++) {
    //for (xx=0;xx<sz;xx++) dataout[sz*sz*ii+sz*yy+xx]=data[szinX*szinY*z+szinX*(t+yy)+(l+xx)];
    memcpy(dataout+(sz*sz*ii+sz*yy+0),data+(szinX*szinY*z+szinX*(t+yy)+(l+0)),sz*sizeof(float));
  }

  leftcoord[ii]=(float) l;
  topcoord[ii]=(float) t;

  return;
}
