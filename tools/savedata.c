/*
 * TODO: implement error-checking
 *
 */

#include "savedata.h"

void saveLocsBinary(char *outFile, int locs) {
  FILE *file = fopen( outFile, "wb" );
  fwrite( (void*)&locs, sizeof(int), 1, file );
  fclose( file ); 
}

void saveArrayBinary( char* outFile, float* data, int size ) {
  FILE *file = fopen( outFile, "wb" );
  fwrite( data, sizeof( float ), size, file );
  fclose( file );
}

void save2DBinary( char* outFile, float* data, int size ) {
  FILE *file = fopen( outFile, "wb" );
  fwrite( data, sizeof( float ), size*size, file );
  fclose( file );
}

void saveArray( char* outFile, float* data, int size ) {
  FILE *file = fopen( outFile, "w" );

  for ( int i = 0; i < size; i ++ ) {
    fprintf( file, "%0.4lf\n", data[i] );
    //fprintf( file, "%d %0.4lf\n", i, data[i] );
  }

  fclose( file );
}

void save2D( char* outFile, float* data1, float* data2, int size ) {
  FILE *file = fopen( outFile, "w" );

  for ( int i = 0; i < size; i ++ ) {
    fprintf( file, "%0.4lf %0.4lf\n", data1[i], data2[i] );
    //fprintf( file, "%d %0.4lf\n", i, data[i] );
  }

  fclose( file );
}

