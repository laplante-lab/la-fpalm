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

#include <memory> // make_unique

class SRscmos {

private:
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

//  std::unique_ptr<float []> imtot;
//  std::unique_ptr<float []> xtot;
//  std::unique_ptr<float []> ytot;
//  std::unique_ptr<float []> ztot;
//  std::unique_ptr<float []> bgtot;
//  std::unique_ptr<float []> lltot;
//  std::unique_ptr<float []> photot;
//  std::unique_ptr<float []> crlbxtot;
//  std::unique_ptr<float []> crlbytot;
//
//  std::unique_ptr<double []> maxprojim;
//  std::unique_ptr<double []> sumim;

  uint32_t total_size;
  uint32_t imsz;
  
  double zm;
  uint32_t imszzm;

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

  uint32_t input_data_width;

  SRscmos();
  SRscmos(const SRscmos &s);
  ~SRscmos();

  // sets the fit parameters
  void setFitParameters(uint32_t min_photon, float PSFSigma, uint32_t zm_color, float SRImage_pixelsize);

  // sets the threshold values
  void setThresholds( uint32_t ll_threshold, uint32_t subregion_cropx, uint32_t subregion_cropy, uint32_t N_photon_min, uint32_t N_photon_max, float loc_uncer_min, float loc_uncer_max);

  // runs FPALM analysis on a given dataset
  void run(NoiseMap* nm, Image<uint16_t> *frame_stack, uint32_t num_dataset);

  // saves reconstruction image to output folder 
  void saveData(const char* output_file_path);

  // maps data between two SRscmos objects
  void mapData(const SRscmos &local_data, uint32_t num_dataset, uint32_t frames_per_worker_fixed, uint32_t remaining_frames);

  void mergeData(float* cp_imtot, float* cp_xtot, float* cp_ytot, float* cp_ztot, float* cp_bgtot, float* cp_lltot, float* cp_photot, float* cp_crlbxtot, float* cp_crlbytot, double* cp_maxprojim, double* cp_sumim, uint32_t cp_total_size, uint32_t num_dataset, uint32_t total_frames, uint32_t remaining_frames, uint32_t data_width);

  // saves data to Matlab data file
  const std::string &inversegain_path = "./noise_map/la-gainmap.mat";

  void writeMatFile(const std::string &gain_path, const std::string &var_path, const std::string &input_path, const std::string &output_path, uint32_t total_frames);

//  // useful for debugging an instance of SRscmos
//  friend std::ostream& operator<<(std::ostream &out, const SRscmos& sr) {
//    out << "total_size: " << sr.total_size;
//    //out << std::setfill(' ') << std::dec << std::setw(7) << cl.index;
//
//    //out << std::hex << setfill(' ') << std::setw(4+(int)(ceil((double)cl.tag_bits/4))) << showbase << cl.tag;
//    out << std::endl;
//    return out;
//  }
};

#endif
