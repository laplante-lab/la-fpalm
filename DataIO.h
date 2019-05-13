#ifndef DATA_IO_H
#define DATA_IO_H

#include <tiffio.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

//#include <inttypes.h>
//#include <math.h>

#include <vector>
#include <algorithm>
#include "dirent.h"

//#include "parameters.h"

//#include <string.h>
//#include <cstdint>

#include "Image.h"

#include <Nd2ReadSdk.h>

#include <experimental/filesystem> // or #include <filesystem>
#include <algorithm>    // std::sort
#include <sys/stat.h>

#include <string>
#include "dirent.h"

// https://github.com/nlohmann/json
#include <nlohmann/json.hpp>

using json = nlohmann::json;

typedef enum{UNKNOWN, ND2, TIF} file_type_t;
//static const char *file_type_names[3] = {"UKNOWN", "ND2", "TIF"};

int createOutput(const std::string &output_path);

void writeImage(const char* filename, double* r, double* g, double* b, uint32_t width, uint32_t height);

string simplifyString(string text);

class DataIO {
private:
  // nd2 metadata
  LIMFILEHANDLE handle;
  json jAttributes;

  std::vector<std::string> image_names; // file_names

  bool is_file(const char *path);
  bool is_dir(const char *path);

  // read_directory()
  // Return an ASCII-sorted vector of filename entries in a given directory.
  // If no path is specified, the current working directory is used.
  //
  // Always check the value of the global 'errno' variable after using this
  // function to see if anything went wrong. (It will be zero if all is well.)
  std::vector <std::string> read_directory(const std::string& path);

public:
  std::string data_path;
  std::string file_name_noext;

  file_type_t input_format;

  uint32_t number_of_frames;
  uint32_t frame_height;
  uint32_t frame_width;

  DataIO(const std::string &data_path);
  ~DataIO();

  std::string getFileName();
  std::string getFileExtension();

  int getNumChannels();
  std::string getChannelName(uint32_t channel_index);

  Image<uint16_t>* parseImageData(uint32_t start_frame, uint32_t end_frame);
  Image<uint16_t>* parseImageData(uint32_t start_frame, uint32_t end_frame, uint32_t channel);
};

#endif

