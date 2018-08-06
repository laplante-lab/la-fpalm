/**
*  @file    Image.h
*  @author  John Ravi (jjravi)
*  @date    6/13/2018
*  @version 2.2
*
*  @brief 2D/3D image stack data structure
*
*  @section DESCRIPTION
*
*  Easy to use data structure to store pixel values. Type is detemined at compile time (see c++ template).
*
*  Usage:
*    Image<type> object_name = new Image<type>(#frames, height, width);
*    pixel = object_name[std::make_tuple(frame #, height #, width #)]
*    pixel = object_name[{height #, width #}]
*    image = object_name[frame #]
*
*/

#ifndef IMAGE_H
#define IMAGE_H

#include <cstdlib>
#include <utility> //std::pair
#include <tuple> // std::tuple

#include <iostream>
#include <iomanip> // std::setw

using namespace std;

template <class E> // Type
class Image {
private:

public:
	E * pixels;
	uint32_t num_frames;
	uint32_t height;
	uint32_t width;

  Image(const uint32_t F, const uint32_t H, const uint32_t W); // constructor
  Image(const uint32_t H, const uint32_t W); // constructor

	~Image(); // deconstructor

  // resize data structure
	void resize(uint32_t F, uint32_t H, uint32_t W);

  // fill all data with zeros
  void fillZeros();

  // delete allocated space
  void freeData();

  // index a pixel using pair, overrides the [] operator
	E& operator [] (const std::pair<uint32_t, uint32_t>& Index) { // get array item 2D array
		return pixels[Index.first * width + Index.second];
		//return pixels[Index.second * width + Index.first];
	}

  // index a pixel using tuple, overrides the [] operator
	E& operator [] (const std::tuple<uint32_t, uint32_t, uint32_t>& Index) { // get array item 3D array
		return pixels[(std::get<0>(Index)*width*height) + (std::get<1>(Index)*width) + std::get<2>(Index)];
	}

  // toString method, overrides the << operator. Use with std::cout
  friend std::ostream& operator<<(std::ostream &out, const Image& img) {

    //out << std::setfill(' ') << std::dec << std::setw(7) << cl.index;
    out << img.num_frames << 'x' << img.height << 'x' << img.width << std::endl;

    for ( uint32_t i = 0; i < img.num_frames; i++ ) {
      out << std::endl << "Frame #" << i << std::endl;
      for ( uint32_t j = 0; j < img.height; j++ ) {
        for ( uint32_t k = 0; k < img.width; k++ ) {
          out << std::setfill(' ') << std::dec << std::setw(3) << (*((Image *)&img))[std::make_tuple(i,j,k)] << " ";
        }
        out << std::endl;
      }
    }
    out << std::endl;
    return out;
  }

};

//int to1D(int x, int y, int z, int xMax, int yMax) {
//  return (z * xMax * yMax) + (y * xMax) + x;
//}

template <class E>
Image<E>::Image(const uint32_t F, const uint32_t H, const uint32_t W) {
	height = H;
	width = W;
	num_frames = F;

	pixels = (E *)malloc(num_frames*height*width * sizeof(E));
}

template <class E>
Image<E>::Image(const uint32_t H, const uint32_t W) {
	height = H;
	width = W;
	num_frames = 1;

	pixels = (E *)malloc(num_frames*height*width * sizeof(E));

}

template <class E>
Image<E>::~Image() {
	if (pixels) {
    free(pixels); // freeing memory
		pixels = NULL;
	}
}

template <class E>
void Image<E>::resize(const uint32_t F, const uint32_t H, const uint32_t W) {

	this->num_frames = F;
	this->height = H;
	this->width = W;

	pixels = (E *)realloc(pixels, num_frames*height*width * sizeof(E));
}

template <class E>
void Image<E>::fillZeros() {
  std::tuple<uint32_t, uint32_t, uint32_t> index;
  for ( uint32_t i = 0; i < num_frames; i++ ) {
    std::get<0>(index)=i;
    for ( uint32_t j = 0; j < height; j++ ) {
      std::get<1>(index)=j;
      for ( uint32_t k = 0; k < width; k++ ) {
        std::get<2>(index)=k;
        pixels[(std::get<0>(index)*width*height) + (std::get<1>(index)*width) + std::get<2>(index)] = 0;
      }
    }
  }
}

template <class E>
void Image<E>::freeData() {
  free(pixels);
}

#endif
