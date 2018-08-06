#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

#include "srhist_color.h"

#include "cHistRecon.h"

#ifdef SERIAL
// NOTHING
#elif OMP
#include <omp.h>
#endif

double* imstretch_linear(int image_size, double* im, int low_in, int high_in, int low_out, int high_out) {

  double* imout = (double*) malloc( image_size*image_size*sizeof(double) );

  int i;

#ifdef OMP
#pragma omp parallel for default(shared)
#endif
  for ( i = 0; i < image_size*image_size; i++ ) {
    imout[i] = (im[i] - low_in) / (high_in - low_in);

    if ( imout[i] > 1 ) {
      imout[i] = 1;
    }
    else if ( imout[i] < 0 ) {
      imout[i] = 0;
    }

    imout[i] = (imout[i]) * (high_out - low_out) + low_out;

  }

  return imout;
}


uint16_t* SRreconstructhist( int sz, int zm, int size, float* xtmp, float* ytmp ) {
  int maxblob=250000;
  int maxk = ceil( (double)size / maxblob ); 

  uint16_t* imtemp = (uint16_t*) calloc( sz*zm * sz*zm, sizeof(uint16_t) );
  
  int bst, bed;

  //printf("maxk: %d\n", maxk);
  int imszzm = sz * zm;

  int szz, index; 
  float* xresult2 = (float*) malloc( sizeof(float) );
  float* yresult2 = (float*) malloc( sizeof(float) );

  int i, j;

  for( int ii = 1; ii <= maxk; ii++ ) { 
    bst = (ii - 1) * maxblob + 1;
    
    if ( ii == maxk ) {
      bed = size;
    }
    else {
      bed = ii * maxblob;
    }

    szz = bed - bst + 1;

    xresult2 = (float*) realloc( xresult2, szz * sizeof(float) );
    yresult2 = (float*) realloc( yresult2, szz * sizeof(float) );

    index = 0;
    for( i = bst-1; i < bed; i++ ) {
      xresult2[index] = xtmp[i] * zm;
      yresult2[index++] = ytmp[i] * zm;
    }

    uint16_t* srim = cHistRecon( imszzm, imszzm, szz, xresult2, yresult2, 0 );

    free( xresult2 );
    free( yresult2 );

    for ( i = 0; i < imszzm; i++ ) {
      for( j = 0; j < imszzm; j++ ) {
        imtemp[i*imszzm + j] = imtemp[i*imszzm + j] + srim[i*imszzm + j];
      }
    }

    free(srim);
  }

  return imtemp;
} 

void srhist_color(int sz, int zm, int size, float* xtot, float* ytot, float* ttot, int segnum, double* rch, double* gch, double* bch) {

  // color map ///
  FILE * color_map_file = fopen("tools/colormapjet.bin","rb");

  if ( !color_map_file ) {
    printf("can not open color_map\n");
    exit(1);
  }

  double *color_map_jet = (double *)malloc( 64 * 3 * sizeof(double) ); // 3 channels, 64
  fread( color_map_jet, sizeof( double ), 64 * 3, color_map_file ); // 64 x 3 double
  fclose (color_map_file);

  double *a = color_map_jet;

  float max = ttot[0]; 
  float min = ttot[0];

  for ( int i = 1; i < size; i++ ) {
    if ( max < ttot[i] )
      max = ttot[i];

    if ( min > ttot[i] )
      min = ttot[i];
  }

  float incre = floor( (max-min) / segnum );
  //printf("incre: %0.4lf\n", incre);
  
  float tst, ted;

  bool* mask;

  int im_size = sz*zm;
  //printf("size: %d\n", im_size);

  float* xtmp;// = (float*) malloc( sizeof(float) );
  float* ytmp;// = (float*) malloc( sizeof(float) );

  int i, j;
  int ii;

  int tmp_size;
  int index;

  uint16_t* tmpim;

#ifdef OMP
#pragma omp parallel for default(shared) private( mask, tst, ted, i, j, tmp_size, index, tmpim, xtmp, ytmp )
#endif
  for( ii = 1; ii <= segnum; ii++ ) {
    tst = (ii - 1)*incre + min;
 
    if( ii == segnum ) {
      ted = max;
    }
    else {
      ted = ii*incre+min;
    }
    
    //printf("tst: %0.4lf ted: %0.4lf\n", tst, ted);

    tmp_size = 0;

    mask = (bool*)malloc( size * sizeof( bool ) );

    for ( i = 0; i < size; i++ ) {
      mask[i] = (ttot[i] >= tst && ttot[i] <= ted);
      if(mask[i]) {
        tmp_size++;
      }
    }

    xtmp = (float*) malloc( tmp_size*sizeof(float) ); 
    ytmp = (float*) malloc( tmp_size*sizeof(float) );

    // xtmp = (float*) realloc( xtmp, tmp_size*sizeof(float) );
    // ytmp = (float*) realloc( ytmp, tmp_size*sizeof(float) );

    index = 0;
    for ( i = 0; i < size; i++ ) {
      if ( mask[i] ) {
        xtmp[index] = xtot[i];
        ytmp[index++] = ytot[i];

        //printf("xtmp: %0.4lf ytmp: %0.4lf\n", xtmp[index-1], ytmp[index-1]);
      }

      if ( index >= tmp_size ) break; // optimization 
    }

    free(mask);

    tmpim = SRreconstructhist(sz,zm,tmp_size,xtmp,ytmp);
    free( xtmp );
    free( ytmp );

    for( i = 0; i < im_size; i++ ) {
      for( j = 0; j < im_size; j++ ) {
#ifdef OMP
#pragma omp atomic
#endif
        rch[i*im_size + j] = rch[i*im_size + j] + ( a[(ii-1)*3 + 0] * (double)tmpim[i*im_size + j] );
#ifdef OMP
#pragma omp atomic
#endif
        gch[i*im_size + j] = gch[i*im_size + j] + ( a[(ii-1)*3 + 1] * (double)tmpim[i*im_size + j] );
#ifdef OMP
#pragma omp atomic
#endif
        bch[i*im_size + j] = bch[i*im_size + j] + ( a[(ii-1)*3 + 2] * (double)tmpim[i*im_size + j] );
      }
    }

    free( tmpim );
  }

  //free( xtmp );
  //free( ytmp );
}
