/**
*  @file    DataIO.cpp
*  @author  John Ravi (jjravi)
*  @date    6/12/2018
*  @version 2.0
*
*  @brief Allow for parsing data from .nd2 files
*
*  @section DESCRIPTION
*
*  TODO: ADD Descriptuon 
*/

#include "DataIO.h"

#include <stdlib.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

//#include <pthread.h>

#include <string>
#include <tuple>

namespace fs = std::experimental::filesystem;
using namespace std;

DataIO::DataIO(const std::string &data_path) {
  this->data_path = data_path;
  input_format = UNKNOWN;

  // Some error checking on input data path
  if ( is_file(data_path.c_str()) ) {
    // check if .nd2 extension
    if(data_path.substr(data_path.find_last_of(".") + 1) == "nd2") {
      input_format = ND2;
      printf("Loading .nd2 file\n");
    }
    else {
      fprintf(stderr, "Input file type is not .nd2\n");
      exit(EXIT_FAILURE);
    }
  }
  else if ( is_dir(data_path.c_str()) ) {
    // data_path is folder, assume TIF
    input_format = TIF;
    printf("Loading .tif images\n");
  }
  else {
    input_format = UNKNOWN;
    fprintf(stderr, "Input data is not found\n");
    exit(EXIT_FAILURE);
  }

  // data loaded here
  switch (input_format) {
    case ND2:
      handle = Lim_FileOpenForReadUtf8(data_path.c_str());

      { // to ensure variable scope for each case
        LIMSTR aAttributes = Lim_FileGetAttributes(handle);
        jAttributes = json::parse(aAttributes);
        Lim_FileFreeString(aAttributes);
      }

      number_of_frames = Lim_FileGetSeqCount(handle);
      frame_height = jAttributes["heightPx"];
      frame_width = jAttributes["widthPx"];
      break;
    case TIF:
      break;
    default:
      break;
  }
}

DataIO::~DataIO() {
  switch (input_format) {
    case ND2:
      if ( handle ) {
        Lim_FileClose(handle);
      }
      break;
    case TIF:
      break;
    default:
      break;
  }
}

std::vector <std::string> DataIO::read_directory(const std::string& path = std::string()) {
	std::vector <std::string> result;
	dirent* de;
	DIR* dp;
	errno = 0;
	dp = opendir(path.empty() ? "." : path.c_str());
	char ext[5];

	if (dp) {
		while (true) {
			errno = 0;
			de = readdir(dp);
			if (de == NULL) break;

			size_t len = strlen(de->d_name);
			strncpy(ext, de->d_name + (len - 4), 5);

			if (strcmp(ext, ".tif") == 0)
				result.push_back(std::string(de->d_name));
		}
		closedir(dp);
		std::sort(result.begin(), result.end());
	}
	return result;
}

int createOutput(const std::string &output_path) {
  // Check if output directory exists
  if (!fs::exists(output_path)) { // Check if src folder exists
    fs::create_directory(output_path); // create src folder 
  }

  return 0;
}

bool DataIO::is_file(const char *path) {
  struct stat buf;
  stat(path, &buf);
  return S_ISREG(buf.st_mode);
}

bool DataIO::is_dir(const char *path) {
  struct stat buf;
  stat(path, &buf);
  return S_ISDIR(buf.st_mode);
}

Image<uint16_t>* DataIO::parseImageData(uint32_t start_frame, uint32_t end_frame) {

  uint32_t subset_size = end_frame - start_frame;
  Image<uint16_t> *frame_stack = new Image<uint16_t>( subset_size, this->frame_height, this->frame_width ); 

  switch (input_format) {
    case ND2:
      {
        LIMPICTURE pPicture;
        Lim_InitPicture( &pPicture, jAttributes["widthPx"], jAttributes["heightPx"], jAttributes["bitsPerComponentSignificant"], jAttributes["bitsPerComponentInMemory"] ); 

        for (uint32_t i = start_frame, ii = 0; i < end_frame; i++, ii++) {
          // get single image data at index i from nd2 file
          Lim_FileGetImageData(handle, i, &pPicture);

          // convert to uint16 2d array
          //uint16_t *pixels = (uint16_t*) pPicture.pImageData;

          //for( uint32_t j = 0; j < this->frame_height; j++ )
          //  for( uint32_t k = 0; k < this->frame_width; k++ )
          //    (*frame_stack)[std::make_tuple(ii, j, k)] = pixels[j*pPicture.uiWidth + k];

          // this line of code is equivalent to the commented nested for loop above
          memcpy( frame_stack->pixels+ii*frame_stack->width*frame_stack->height, pPicture.pImageData, this->frame_height*this->frame_width*sizeof(uint16_t) ); 
        }
        Lim_DestroyPicture(&pPicture);
      }
      break;
    case TIF:
      break;
    default:
      break;
  }

  return frame_stack;
}


/*
 * This routine writes an RGB TIFF image to "filename" scanline by scanline.
 * The arguments to the function are the filename, the raster data in ABGR
 * format (8 bits per color), the number of rows, and the number of columns.
 * This function automatically builds the TIFF information for a single RGB
 * image of the given size.
 */

#define INTEL 1

static int tiff_raster2file(char filename[], uint32 *raster, uint32 nrow, uint32 ncol) {
  uint16_t i, j;
  long stripsperimage = 0;
  long rasterbase;
  unsigned char *tmpStorage;
  unsigned char *tmpPtr;
  TIFF* tif = TIFFOpen(filename, "w");
  if (tif) {
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, ncol);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, nrow);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, 1); /* none */
    TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
    TIFFSetField(tif, TIFFTAG_XRESOLUTION, 1200.0);
    TIFFSetField(tif, TIFFTAG_YRESOLUTION, 1200.0);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

    /* calculate the number of strips we want in the TIFF file */
    stripsperimage = nrow * ncol / 8000;
    if (stripsperimage < 1)
      stripsperimage = 1;
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, stripsperimage);

    tmpStorage = (unsigned char *)malloc(sizeof(unsigned char) * ncol * 3);
    if (!tmpStorage) {
      fprintf(stderr, "Memory allocation error\n");
      TIFFClose(tif);
      return(-1);
    }
    if (raster != NULL) {
      for (i = 0; i<nrow; i++) {
        rasterbase = (nrow - i - 1) * ncol;
        for (j = 0; j<ncol * 3;) {
          tmpPtr = (unsigned char *)&(raster[rasterbase]);
#if INTEL
          tmpStorage[j] = tmpPtr[0];
          j++;
          tmpStorage[j] = tmpPtr[1];
          j++;
          tmpStorage[j] = tmpPtr[2];
          j++;
#else
          tmpStorage[j] = tmpPtr[3];
          j++;
          tmpStorage[j] = tmpPtr[2];
          j++;
          tmpStorage[j] = tmpPtr[1];
          j++;
#endif
          rasterbase++;
        }
        if (TIFFWriteScanline(tif, tmpStorage, i, 0) != 1) {
          free(tmpStorage);
          fprintf(stderr, "Unable to write scanline %d\n", i);
          return(-1);
        }
      }

      free(tmpStorage);
      TIFFClose(tif);
      return(0);
    }
    else {
      TIFFClose(tif);
      free(tmpStorage);
      fprintf(stderr, "image is NULL\n");
      return(-1);
    }
  }
  else {
    fprintf(stderr, "Unable to open %s for writing\n", filename);
    return(-1);
  }
}

//static void transpose(double *mat, int size) {
//  for (int n = 0; n < (size - 2); n++)
//    for (int m = n + 1; m < (size - 1); m++)
//      mat[n*size + m] = mat[m*size + n];
//}

// An Inplace function to rotate a N x N matrix 
// by 90 degrees in anti-clockwise direction 
static void rotateMatrix(double *mat, int N) {
  // Consider all squares one by one     
  for (int x = 0; x < N / 2; x++) {
    // Consider elements in group of 4 in         
    // current square         
    for (int y = x; y < N - x - 1; y++) {
      // store current cell in temp variable
      int temp = mat[x*N + y];
      // move values from right to top
      mat[x*N + y] = mat[y*N + (N - 1 - x)];
      // move values from bottom to right
      mat[y*N + (N - 1 - x)] = mat[(N - 1 - x)*N + (N - 1 - y)];
      // move values from left to bottom
      mat[(N - 1 - x)*N + (N - 1 - y)] = mat[(N - 1 - y)*N + x];
      // assign temp to left
      mat[(N - 1 - y)*N + x] = temp;
    }
  }
}

void writeImage(const char* filename, double* r, double* g, double* b, uint32_t width, uint32_t height) {

  uint32_t* raster = (uint32_t*)malloc(width*height * sizeof(uint32_t));
  unsigned char* color = (unsigned char*)malloc(3 * sizeof(unsigned char));

  rotateMatrix(r, height);
  rotateMatrix(g, height);
  rotateMatrix(b, height);

  for (uint32_t i = 0; i < width; i++) {
    for (uint32_t j = 0; j < height; j++) {

      color[0] = (unsigned char)round(r[i*width + j]);
      color[1] = (unsigned char)round(g[i*width + j]);
      color[2] = (unsigned char)round(b[i*width + j]);

      unsigned char* one = (unsigned char *)&raster[i*width + j];
      one[0] = color[0];
      one[1] = color[1];
      one[2] = color[2];
    }
  }

  tiff_raster2file((char*)filename, raster, width, height);
  free(raster);
  free(color);
}

