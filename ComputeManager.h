#ifndef COMPUTE_MANAGER_H
#define COMPUTE_MANAGER_H

#include "SRscmos.h"    // FPALM analysis kernels
#include "DataIO.h"     // data parser (.nd2 support and .tiff support)
#include "tools/plot.h" // data plot

#ifdef OMP
  #include <omp.h>      // cpu thread level parallelism
#endif

#include <stdbool.h>
#include <math.h> // ceil
#include <string.h> // to_string 

#include "Image.h"

class ComputeManager {

private:
  DataIO *input;

public:

  std::string inversegain_path;
  std::string varmap_path;

  //// default parameters ////
  // Fit Parameters
  uint32_t min_photon;
  float PSFSigma;
  uint32_t zm_color;
  float SRImage_pixelsize;

  // Thresholding
  uint32_t ll_threshold;
  uint32_t subregion_cropx;
  uint32_t subregion_cropy;
  uint32_t N_photon_min;
  uint32_t N_photon_max;
  float loc_uncer_min;
  float loc_uncer_max;
  ////////////////////////////

  ComputeManager();
  ~ComputeManager();

  void compute(const std::string &data_path, const std::string &output_path, int new_start=-1, int new_end=-1);

  //void compute(const std::string &data_path, const std::string &output_path);
  //void compute(const std::string inversegain_path, const std::string varmap_path, const std::string data_path, const std::string output_path, 
  //              int min_photon, float PSFSigma, int zm_color, float SRImage_pixelsize, 
  //              int ll_threshold, int subregion_cropx, int subregion_cropy, int N_photon_min, int N_photon_max, float loc_uncer_min, float loc_uncer_max,
  //              int subset_size);

  //void compute2(const std::string inversegain_path, const std::string varmap_path, const std::string data_path, const std::string output_path, 
  //               int min_photon, float PSFSigma, int zm_color, float SRImage_pixelsize, 
  //              int ll_threshold, int subregion_cropx, int subregion_cropy, int N_photon_min, int N_photon_max, float loc_uncer_min, float loc_uncer_max,
  //               int new_start, int new_end);

//  SRscmos* compute_subset(char* inversegain_path, char* varmap_path, char* data_path, char* output_path, 
//                          int min_photon, float PSFSigma, int zm_color, float SRImage_pixelsize, 
//                          int ll_threshold, int subregion_cropx, int subregion_cropy, int N_photon_min, int N_photon_max, float loc_uncer_min, float loc_uncer_max,
//                          int new_start, int new_end);
};

#endif

