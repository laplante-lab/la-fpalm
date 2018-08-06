#ifndef SR_SCMOS_H
#define SR_SCMOS_H

#include "parameters.h"
#include "noise_map/NoiseMap.h"
#include "gaussmle/kernelCall.h" 

extern "C" {
#include "tools/cMakeSubregions.h"
#include "tools/cHistRecon.h"

#include "tools/srhist_color.h"

#include "tools/savedata.h"
}

#include "tools/convolution.h"

#include "DataIO.h"
#include <math.h>

#include "Image.h"

typedef double DoubleImage[253][253];

class SRscmos {

private:
  void findRegionsOfInterest(DoubleImage* unifim, int total_frames, bool(*im_max)[253][253]);
  Image<double>* varunif(Image<double> *tmp1, NoiseMap* nm);

public:
  float* imtot;
  float* xtot;
  float* ytot;
  float* ztot;
  float* bgtot;
  float* lltot;
  float* photot;
  float* crlbxtot;
  float* crlbytot;

  double* maxprojim;
  double* sumim;

  int total_size;
  int imsz;
  
  double zm;
  int imszzm;

  // Fit Parameters
  int min_photon;
  float PSFSigma;
  int zm_color;
  float SRImage_pixelsize;

  // Thresholding
  int ll_threshold;
  int subregion_cropx;
  int subregion_cropy;
  int N_photon_min;
  int N_photon_max;
  float loc_uncer_min;
  float loc_uncer_max;

  uint32_t input_data_width;

  SRscmos();
  ~SRscmos();

  // sets the fit parameters
  void setFitParameters(int min_photon, float PSFSigma, int zm_color, float SRImage_pixelsize);

  // sets the threshold values
  void setThresholds( int ll_threshold, int subregion_cropx, int subregion_cropy, int N_photon_min, int N_photon_max, float loc_uncer_min, float loc_uncer_max);

  // runs FPALM analysis on a given dataset
  void run(NoiseMap* nm, Image<uint16_t> *frame_stack, int num_dataset);

  // runs FPALM analysis on a given dataset
//  void run(NoiseMap* nm, Frame256x256* images, int num_dataset, int total_frames);

  // saves reconstruction image to output folder 
  void saveData(const char* output_file_path);

  // merges data between two SRscmos objects
  void mergeData(float* cp_imtot, float* cp_xtot, float* cp_ytot, float* cp_ztot, float* cp_bgtot, float* cp_lltot, float* cp_photot, float* cp_crlbxtot, float* cp_crlbytot, double* cp_maxprojim, double* cp_sumim, int cp_total_size, int num_dataset, int total_frames, int remaining_frames);

  // merges data between two SRscmos objects
  void mergeData(float* imtot, float* xtot, float* ytot, float* ztot, float* bgtot, float* lltot, float* photot, float* crlbxtot, float* crlbytot, double* maxprojim, double* sumim, int cp_total_size, int num_dataset, int total_frames, int remaining_frames, int number_of_iterations);

  // saves data to Matlab data file
  const std::string &inversegain_path = "./noise_map/la-gainmap.mat";

  void writeMatFile(const std::string &gain_path, const std::string &var_path, const std::string &input_path, const std::string &output_path, int total_frames);

};

#endif
