#ifndef SAVE_DATA_H
#define SAVE_DATA_H

#include <stdio.h>

void saveLocsBinary(char *outFile, int locs);
void saveArrayBinary( char* outFile, float* data, int size );
void save2DBinary( char* outFile, float* data, int size );

void saveArray( char* outFile, float* data, int size );
void save2D( char* outFile, float* data1, float* data2, int size );

#endif
