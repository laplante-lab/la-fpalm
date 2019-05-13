#include "diplib.h"

#include "diplib/linear.h"
#include "diplib/morphology.h"
#include "diplib/multithreading.h"

#include "SRscmos.h"

//#include "omp.h"

int to1D(int x, int y, int z, int xMax, int yMax) {
  return (z * xMax * yMax) + (y * xMax) + x;
}

SRscmos::SRscmos() {
  input_data_width = 0;

  // Fit Parameters
  min_photon = 200;
  PSFSigma = 1.3;
  zm_color = 20;
  SRImage_pixelsize = 0.1;

  // Thresholding
  ll_threshold = 80;
  subregion_cropx = 1;
  subregion_cropy = 1;
  N_photon_min = 60;
  N_photon_max = 6000;
  loc_uncer_min = 0.01;
  loc_uncer_max = 2;

  total_size = 0;

  xtot = (float *)malloc(sizeof(float));
  ytot = (float *)malloc(sizeof(float));
  ztot = (float *)malloc(sizeof(float));
  bgtot = (float *)malloc(sizeof(float));
  lltot = (float *)malloc(sizeof(float));
  photot = (float *)malloc(sizeof(float));
  crlbxtot = (float *)malloc(sizeof(float));
  crlbytot = (float *)malloc(sizeof(float));
}

SRscmos::SRscmos(const SRscmos &s) {
  this->input_data_width = s.input_data_width;

  // Fit Parameters
  this->min_photon = s.min_photon;
  this->PSFSigma = s.PSFSigma;
  this->zm_color = s.zm_color;
  this->SRImage_pixelsize = s.SRImage_pixelsize;

  // Thresholding
  this->ll_threshold = s.ll_threshold;
  this->subregion_cropx = s.subregion_cropx;
  this->subregion_cropy = s.subregion_cropy;
  this->N_photon_min = s.N_photon_min;
  this->N_photon_max = s.N_photon_max;
  this->loc_uncer_min = s.loc_uncer_min;
  this->loc_uncer_max = s.loc_uncer_max;

  this->imsz = s.imsz;
  this->zm = s.zm;
  this->imszzm = s.imszzm;

  imtot = (float *)calloc(imszzm * imszzm, sizeof(float));
  xtot = (float *)malloc(sizeof(float));
  ytot = (float *)malloc(sizeof(float));
  ztot = (float *)malloc(sizeof(float));
  bgtot = (float *)malloc(sizeof(float));
  lltot = (float *)malloc(sizeof(float));
  photot = (float *)malloc(sizeof(float));
  crlbxtot = (float *)malloc(sizeof(float));
  crlbytot = (float *)malloc(sizeof(float));

  maxprojim = (double *)calloc(imsz * imsz, sizeof(double));
  sumim = (double *)calloc(imsz * imsz, sizeof(double));

//  this->imtot = std::make_unique<float []>(s.imszzm * s.imszzm);
//  this->maxprojim = std::make_unique<double[]>(s.imsz*s.imsz);
//  this->sumim = std::make_unique<double[]>(s.imsz*s.imsz);
//
//  xtot = std::make_unique<float []>(s.total_size);
//  ytot = std::make_unique<float []>(s.total_size);
//  ztot = std::make_unique<float []>(s.total_size);
//
//  bgtot = std::make_unique<float []>(s.total_size);
//  lltot = std::make_unique<float []>(s.total_size);
//  photot = std::make_unique<float []>(s.total_size);
//  crlbxtot = std::make_unique<float []>(s.total_size);
//  crlbytot = std::make_unique<float []>(s.total_size);

  this->total_size = 0;
}



void SRscmos::setFitParameters(uint32_t min_photon, float PSFSigma, uint32_t zm_color, float SRImage_pixelsize) {
  this->min_photon = min_photon;
  this->PSFSigma = PSFSigma;
  this->zm_color = zm_color;
  this->SRImage_pixelsize = SRImage_pixelsize;
}

void SRscmos::setThresholds( uint32_t ll_threshold, uint32_t subregion_cropx, uint32_t subregion_cropy, uint32_t N_photon_min, uint32_t N_photon_max, float loc_uncer_min, float loc_uncer_max ) {
  this->ll_threshold = ll_threshold;
  this->subregion_cropx = subregion_cropx;
  this->subregion_cropy = subregion_cropy;
  this->N_photon_min = N_photon_min;
  this->N_photon_max = N_photon_max;
  this->loc_uncer_min = loc_uncer_min;
  this->loc_uncer_max = loc_uncer_max;
}

SRscmos::~SRscmos() {
  free(xtot);
  free(ytot);
  free(ztot);
  free(bgtot);
  free(lltot);
  free(photot);
  free(crlbxtot);
  free(crlbytot);

  free(imtot);

  free(maxprojim);
  free(sumim);
}

void SRscmos::run(NoiseMap* nm, Image<uint16_t> *frame_stack, uint32_t num_dataset) {
  this->input_data_width = frame_stack->width;

  ///////////////
  imsz = frame_stack->width - 3;

  zm = 1 / (double)SRImage_pixelsize;
  imszzm = (uint32_t)(imsz*zm);

  imtot = (float *)calloc(imszzm * imszzm, sizeof(float));
  maxprojim = (double *)calloc(imsz * imsz, sizeof(double));
  sumim = (double *)calloc(imsz * imsz, sizeof(double));

  //imtot = std::make_unique<float []>(imszzm * imszzm);
  //maxprojim = std::make_unique<double[]>(imsz*imsz);
  //sumim = std::make_unique<double[]>(imsz*imsz);
  ///////////////

  // This applies the noise filter to the input images
  Image<double> *fin = new Image<double>( frame_stack->num_frames, frame_stack->height, frame_stack->width ); 

  uint32_t mapbase = 0;

  // crop out subsection of camera noise
  if ( frame_stack->height < nm->gain_map->height ) {
    mapbase = (nm->gain_map->height - frame_stack->height)/2;
  }

  //printf("mapbase: %d\n", mapbase);

  for (uint32_t i = 0; i < frame_stack->num_frames; i++)
    for (uint32_t j = 0; j < frame_stack->height; j++)
      for (uint32_t k = 0; k < frame_stack->width; k++)
        (*fin)[std::make_tuple(i, j, k)] = ((double)(*frame_stack)[std::make_tuple(i, j, k)] - 
                      (*nm->o_map)[{k+mapbase, j+mapbase}])/ (*nm->gain_map)[{k+mapbase,j+mapbase}];


  Image<float> *fcin = new Image<float>( frame_stack->num_frames, frame_stack->height-3, frame_stack->width-3 ); 
  Image<double> *tmp1 = new Image<double>( frame_stack->num_frames, frame_stack->height-3, frame_stack->width-3 ); 

  for (uint32_t i = 0; i < fin->num_frames; i++)
    for (uint32_t j = 0; j < fin->height-3; j++)
      for (uint32_t k = 0; k < fin->width-3; k++)
        (*fcin)[std::make_tuple(i,k,j)] = (float)(*fin)[std::make_tuple(i,j+1,k+1)];

  delete fin;

  for (uint32_t i = 0; i < fcin->num_frames; i++)
    for (uint32_t j = 0; j < fcin->height; j++)
      for (uint32_t k = 0; k < fcin->width; k++)
        (*tmp1)[std::make_tuple(i,j,k)] = (*fcin)[std::make_tuple(i,k,j)]/(*nm->var_map)[{k+mapbase+1,j+mapbase+1}];


  double *kernel = (double *)malloc(9 * sizeof(double));
  for (uint32_t i = 0; i < 9; i++) {
    kernel[i] = 1;
  }

  ///////// 2D convolution kernel size 3x3 /////////////
  Image<double> *cim = new Image<double>(tmp1->num_frames, tmp1->height, tmp1->width); 
  for (uint32_t i = 0; i < tmp1->num_frames; i++) {
    convolve2DSeparable( (tmp1->pixels + i*tmp1->width*tmp1->height), 
                         (cim->pixels + i*cim->width*cim->height), tmp1->height, tmp1->width,
                         kernel, 3, kernel, 3);
  }

  Image<double> *vv = new Image<double>(cim->num_frames, cim->height, cim->width);
  for (uint32_t i = 0; i < vv->num_frames; i++)
    for (uint32_t j = 0; j < vv->height; j++)
      for (uint32_t k = 0; k < vv->width; k++)
        (*vv)[std::make_tuple(i,j,k)] = 1/(*nm->var_map)[{k+mapbase+1,j+mapbase+1}];

  Image<double> *wim = new Image<double>(vv->num_frames, vv->height, vv->width);
  for (uint32_t i = 0; i < wim->num_frames; i++)
    convolve2DSeparable( (vv->pixels + i*vv->width*vv->height), 
                         (wim->pixels + i*wim->width*wim->height), wim->height, wim->width,
                         kernel, 3, kernel, 3);


  Image<double> *filteredim1 = new Image<double>(cim->num_frames, cim->height, cim->width);
  for (uint32_t i = 0; i < filteredim1->num_frames; i++)
    for (uint32_t j = 0; j < filteredim1->height; j++)
      for (uint32_t k = 0; k < filteredim1->width; k++)
        (*filteredim1)[std::make_tuple(i,j,k)] = (*cim)[std::make_tuple(i,j,k)] / (*wim)[std::make_tuple(i,j,k)];

  /////////// 2D convolution kernel size 9x9 ///////////////
  for (uint32_t i = 0; i < tmp1->num_frames; i++) {
    convolve2DSeparable( (tmp1->pixels + i*tmp1->width*tmp1->height), 
                         (cim->pixels + i*cim->width*cim->height), tmp1->height, tmp1->width,
                         kernel, 9, kernel, 9);
  }

  delete tmp1;

  for (uint32_t i = 0; i < vv->num_frames; i++) {
    convolve2DSeparable( (vv->pixels + i*vv->width*vv->height), 
                         (wim->pixels + i*wim->width*wim->height), vv->height, vv->width,
                         kernel, 9, kernel, 9);
  }

  delete vv;

  Image<double> *filteredim2 = new Image<double>(cim->num_frames, cim->height, cim->width);

  for (uint32_t i = 0; i < filteredim2->num_frames; i++)
    for (uint32_t j = 0; j < filteredim2->height; j++)
      for (uint32_t k = 0; k < filteredim2->width; k++)
        (*filteredim2)[std::make_tuple(i,j,k)] = (*cim)[std::make_tuple(i,j,k)] / (*wim)[std::make_tuple(i,j,k)];

  delete cim;
  delete wim;

  Image<double> *unifim = new Image<double>(filteredim1->num_frames, filteredim1->height, filteredim1->width);
  for (uint32_t i = 0; i < unifim->num_frames; i++)
    for (uint32_t j = 0; j < unifim->height; j++)
      for (uint32_t k = 0; k < unifim->width; k++)
        (*unifim)[std::make_tuple(i,j,k)] = (*filteredim1)[std::make_tuple(i,j,k)] - (*filteredim2)[std::make_tuple(i,j,k)];

  free(kernel);

  //printf("convolution done\n");

  for (uint32_t i = 0; i < unifim->num_frames; i++)
    for (uint32_t j = 0; j < unifim->height; j++)
      for (uint32_t k = 0; k < unifim->width; k++)
        (*unifim)[std::make_tuple(i,j,k)] = (*filteredim1)[std::make_tuple(i,j,k)] - (*filteredim2)[std::make_tuple(i,j,k)];

  delete filteredim1;
  delete filteredim2;

  //delete frame_stack;

  Image<bool> *im_max = new Image<bool>(unifim->num_frames, unifim->height, unifim->width);

  double thresh = this->min_photon / 2 / pi / PSFSigma / PSFSigma / 3;
  //printf("thresh: %0.6lf\n", thresh);

  dip::SetNumberOfThreads(1); // diplib has broken parallel mp support

  for (uint32_t ii = 0; ii < unifim->num_frames; ii++) {

    dip::Image img(
      dip::NonOwnedRefToDataSegment((void *)(unifim->pixels+ii*unifim->width*unifim->height)),
      (void *)(unifim->pixels+ii*unifim->width*unifim->height),
      dip::DT_DFLOAT,
      dip::UnsignedArray{ unifim->height, unifim->width },
      { 1, unifim->width }
    );

    dip::Image out;
    dip::Dilation(img, out, { { 5,5 }, "rectangular" }, { "add min" });
    
    double *maxf = (double *)(out.Data());
 
    for (uint32_t i = 0; i < im_max->height; i++)
      for (uint32_t j = 0; j < im_max->width; j++)
        (*im_max)[std::make_tuple(ii,i,j)] = ((*unifim)[std::make_tuple(ii,i,j)] >= 0.999*maxf[i*im_max->width + j]) && ((*unifim)[std::make_tuple(ii,i,j)] > thresh);
  }

  delete unifim;

  //printf("dilation done\n");

  //////////////////
  double *a = (double *)malloc(sizeof(double));
  
  uint32_t num_regions = 0;

  for (uint32_t it = 0; it < im_max->num_frames; it++)
    for (uint32_t j = 0; j < im_max->height; j++)
      for (uint32_t i = 0; i < im_max->width; i++)
        if ((*im_max)[std::make_tuple(it,i,j)]) {
          num_regions++;
          a = (double *)realloc(a, num_regions * sizeof(double));
          a[num_regions - 1] = to1D(i, j, it, im_max->height, im_max->width);
        }

  auto x = std::make_unique<double[]>(num_regions);
  auto y = std::make_unique<double[]>(num_regions);
  auto z = std::make_unique<double[]>(num_regions);

  uint32_t imsz = im_max->width;

  //printf("num_regions: %d\n", num_regions);

  for (uint32_t i = 0; i < num_regions; i++) {
    z[i] = floor(a[i] / im_max->height / im_max->width);
    x[i] = fmod(fmod(a[i], (imsz*imsz)), imsz);
    y[i] = floor(fmod(a[i], (imsz*imsz)) / imsz);
    //printf("%0.4lf  %0.4lf  %0.4lf\n", x[i], y[i], z[i]);
  }

  delete im_max;
  free(a);

  //////////////////
  // sidemask //
  auto sidemask = std::make_unique<bool[]>(num_regions);
  uint32_t new_num_regions = num_regions;

  for (uint32_t i = 0; i < num_regions; i++) {
    sidemask[i] = (x[i] == 0) || (x[i] == frame_stack->width - 2 - 1) || (y[i] == 0) || (y[i] == frame_stack->height - 2 - 1);
    if (sidemask[i]) {
      new_num_regions--;
    }
  }

  auto temp_x = std::make_unique<double[]>(new_num_regions);
  auto temp_y = std::make_unique<double[]>(new_num_regions);
  auto temp_z = std::make_unique<double[]>(new_num_regions);

  uint32_t ind = 0;
  for (uint32_t i = 0; i < num_regions; i++) {
    if (!sidemask[i]) {
      temp_x[ind] = x[i];
      temp_y[ind] = y[i];
      temp_z[ind++] = z[i];
      //printf("%3lf %3lf %3lf\n", x[i], y[i], z[i]);
    }
  }

  x = std::move(temp_x);
  y = std::move(temp_y);
  z = std::move(temp_z);
  num_regions = new_num_regions;
  ////////////////////////////

  ///////// isolate sub-regions /////////////
  //printf("first: x[0]:%0.1lf    y[0]:%0.1lf    z[0]:%0.1lf\n", x[0], y[0], z[0] );
  //printf("center: %0.4lf\n", (*fin)[std::make_tuple((int)z[0],(int)x[0],(int)y[0])] );
  //printf("num_regions:%d\n", num_regions);

  auto t2 = std::make_unique<float[]>(num_regions);
  auto l2 = std::make_unique<float[]>(num_regions);
  
  Image<float> *varmap_cut = new Image<float>(fcin->num_frames, fcin->height, fcin->width);
  Image<float> *gainmap_cut = new Image<float>(fcin->num_frames, fcin->height, fcin->width);

  for (uint32_t i = 0; i < fcin->num_frames; i++)
    for (uint32_t j = 0; j < fcin->height; j++)
      for (uint32_t k = 0; k < fcin->width; k++) {
        (*varmap_cut)[std::make_tuple(i,k,j)] = (float)(*(nm->var_map))[{k+mapbase+1,j+mapbase+1}];
        (*gainmap_cut)[std::make_tuple(i,k,j)] = (float)(*(nm->gain_map))[{k+mapbase+1,j+mapbase+1}];
      }

  float *subregion = cMakeSubregions(x.get(), y.get(), z.get(), num_regions, FIT_BOX_SIZE, fcin->pixels, l2.get(), t2.get(), fcin->height, fcin->width, fcin->num_frames );

  auto tmp_l = std::make_unique<float[]>(num_regions);
  auto tmp_t = std::make_unique<float[]>(num_regions);
 
  float *subvarim = cMakeSubregions(x.get(), y.get(), z.get(), num_regions, FIT_BOX_SIZE, varmap_cut->pixels, tmp_l.get(), tmp_t.get(), fcin->height, fcin->width, fcin->num_frames);
  float *subgainim = cMakeSubregions(x.get(), y.get(), z.get(), num_regions, FIT_BOX_SIZE, gainmap_cut->pixels, tmp_l.get(), tmp_t.get(), fcin->height, fcin->width, fcin->num_frames);

  delete varmap_cut;
  delete gainmap_cut;

  // find the lowres max 
  for (uint32_t i = 0; i < fcin->num_frames; i++)
    for (uint32_t j = 0; j < fcin->height; j++)
      for (uint32_t k = 0; k < fcin->width; k++) {
        if((*fcin)[std::make_tuple(i,k,j)] > maxprojim[j*fcin->width + k]) {
          maxprojim[j*fcin->width + k] = (*fcin)[std::make_tuple(i,k,j)];
        }
      }

  // find the lowres sum
  for (uint32_t i = 0; i < fcin->num_frames; i++)
    for (uint32_t j = 0; j < fcin->height; j++)
      for (uint32_t k = 0; k < fcin->width; k++) {
        sumim[j*fcin->width + k] = sumim[j*fcin->width + k] + (*fcin)[std::make_tuple(i,k,j)];
      }

  delete fcin;

  //////// SCMOS fit ///////
  printf("floats size: %d\n", num_regions);

  auto sub_x = std::make_unique<float[]>(num_regions);
  auto sub_y = std::make_unique<float[]>(num_regions);
  auto sub_photon = std::make_unique<float[]>(num_regions);
  auto sub_background = std::make_unique<float[]>(num_regions);
  auto CRLB2 = std::make_unique<float[]>(4*num_regions);
  auto LL2 = std::make_unique<float[]>(num_regions);

  printf("PSFSigma: %0.4lf, num_dataset: %d \n", PSFSigma, num_dataset);

  SRsCMOS_MLE(num_dataset, num_regions, (float *)subregion, (float *)subvarim, (float *)subgainim, FIT_TYPE, PSFSigma, ITERATIONS, sub_x.get(), sub_y.get(), sub_photon.get(), sub_background.get(), CRLB2.get(), LL2.get());

  free(subregion);
  free(subvarim);
  free(subgainim);

  printf("finished GPU, %d\n", num_dataset);

  auto llmask2 = std::make_unique<bool[]>(num_regions);
  auto xymask = std::make_unique<bool[]>(num_regions);
  auto intmask = std::make_unique<bool[]>(num_regions);

  for (uint32_t i = 0; i < num_regions; i++)
    llmask2[i] = (LL2[i] * (-2)) > ll_threshold;

  for (uint32_t i = 0; i < num_regions; i++)
    xymask[i] = (sub_x[i] >(FIT_BOX_SIZE - subregion_cropx - 1)) || (sub_x[i] < subregion_cropx) || (sub_y[i] > (FIT_BOX_SIZE - subregion_cropx - 1)) || (sub_y[i] < subregion_cropx);

  for (uint32_t i = 0; i < num_regions; i++)
    intmask[i] = (sub_photon[i] < N_photon_min) || (sub_photon[i] > N_photon_max);

  auto uncermask = std::make_unique<bool[]>(num_regions);

  for (uint32_t i = 0; i < num_regions; i++)
    uncermask[i] = (sqrt(CRLB2[i + 1 * num_regions]) < loc_uncer_min) || (sqrt(CRLB2[i + 0 * num_regions]) > loc_uncer_max);

  auto totmask = std::make_unique<bool[]>(num_regions);

  uint32_t mask_size = 0;
  for (uint32_t i = 0; i < num_regions; i++) {
    totmask[i] = llmask2[i] || xymask[i] || intmask[i] || uncermask[i];

    if (!totmask[i])
      mask_size++;
  }

  uint32_t new_size = mask_size;
  printf("new_size: %d\n", new_size);

  auto mask_sub_x = std::make_unique<float[]>(new_size);
  auto mask_sub_y = std::make_unique<float[]>(new_size);
  auto mask_sub_photon = std::make_unique<float[]>(new_size);
  auto mask_sub_background = std::make_unique<float[]>(new_size);
  auto mask_CRLB2 = std::make_unique<float[]>(4*new_size);
  auto mask_l2 = std::make_unique<float[]>(new_size);
  auto mask_t2 = std::make_unique<float[]>(new_size);
  auto mask_z = std::make_unique<double[]>(new_size);
  auto mask_LL2 = std::make_unique<float[]>(new_size);

  uint32_t indexer = 0;
  for (uint32_t i = 0; i < num_regions; i++) {
    if (!totmask[i]) {
      mask_sub_x[indexer] = sub_x[i];
      mask_sub_y[indexer] = sub_y[i];
      mask_sub_photon[indexer] = sub_photon[i];
      mask_sub_background[indexer] = sub_background[i];
      for (uint32_t j = 0; j < 4; j++)
        mask_CRLB2[indexer + j*new_size] = CRLB2[i + j*num_regions];
      mask_l2[indexer] = (float)l2[i];
      mask_t2[indexer] = t2[i];
      mask_z[indexer] = z[i];
      mask_LL2[indexer] = LL2[i];

      indexer++;
    }
  }

  sub_x = std::move(mask_sub_x);
  sub_y = std::move(mask_sub_y);
  sub_photon = std::move(mask_sub_photon);
  sub_background = std::move(mask_sub_background);
  CRLB2 = std::move(mask_CRLB2);
  l2 = std::move(mask_l2);
  t2 = std::move(mask_t2);
  z = std::move(mask_z);
  LL2 = std::move(mask_LL2);
  
  //printf("start reconstruction\n");
  //// Reconstruction /////
  double maxblob = 100000;
  uint32_t maxk = (uint32_t)ceil(new_size / maxblob);
  //printf("maxk:%d\n", maxk);
 
  ////////CONTINUE/////////
  float *xresult2, *yresult2;

  float* xtmp = (float*)malloc(sizeof(float));
  float* ytmp = (float*)malloc(sizeof(float));
  uint32_t xy_tmp_size = 0;

  /////// Parse data /////////////////
  for (uint32_t ii = 0; ii < maxk; ii++) {
    uint32_t bst = (ii)*maxblob + 1;
    //printf("bst: %d\n", bst);

    uint32_t bed;

    if (ii == (maxk - 1))
      bed = new_size;
    else
      bed = ii * maxblob;

    xresult2 = (float *)malloc((bed - bst + 1) * sizeof(float));
    yresult2 = (float *)malloc((bed - bst + 1) * sizeof(float));
 
    uint32_t indx = 0;
    for (uint32_t i = bst - 1; i < bed; i++) {
      xresult2[indx] = sub_y[i] + l2[i];
      yresult2[indx++] = sub_x[i] + t2[i];
    }

    //printf("bed-bst: %d\n", (bed-bst+1));
    //printf("indx: %d\n", indx);
    //printf("imszzm: %d\n", imszzm);

    float* tempx = (float*)malloc((bed - bst + 1) * sizeof(float));
    float* tempy = (float*)malloc((bed - bst + 1) * sizeof(float));
  
    for (uint32_t i = 0; i < indx; i++) {
      tempx[i] = xresult2[i] * zm;
      tempy[i] = yresult2[i] * zm;
    }
  
    xy_tmp_size = xy_tmp_size + indx;
    xtmp = (float*)realloc(xtmp, xy_tmp_size * sizeof(float));
    ytmp = (float*)realloc(ytmp, xy_tmp_size * sizeof(float));
  
    uint16_t *im2 = cHistRecon(imszzm, imszzm, indx, tempx, tempy, 1);
 
    for (uint32_t i = 0; i < imszzm; i++)
      for (uint32_t j = 0; j < imszzm; j++)
        imtot[i*imszzm + j] = imtot[i*imszzm + j] + im2[i*imszzm + j];

    memcpy((void*)(xtmp + xy_tmp_size - indx), (void*)xresult2, (indx) * sizeof(float));
    memcpy((void*)(ytmp + xy_tmp_size - indx), (void*)yresult2, (indx) * sizeof(float));

    free(im2);

    free(tempx);
    free(tempy);

    free(xresult2);
    free(yresult2);
  }
  
  // append results to previous data
  total_size = total_size + xy_tmp_size;

  //TODO: can be optimized
  xtot = (float *)realloc(xtot, total_size * sizeof(float));
  ytot = (float *)realloc(ytot, total_size * sizeof(float));
  ztot = (float *)realloc(ztot, total_size * sizeof(float));

  bgtot = (float *)realloc(bgtot, total_size * sizeof(float));
  lltot = (float *)realloc(lltot, total_size * sizeof(float));
  photot = (float *)realloc(photot, total_size * sizeof(float));
  crlbxtot = (float *)realloc(crlbxtot, total_size * sizeof(float));
  crlbytot = (float *)realloc(crlbytot, total_size * sizeof(float));

//  xtot = std::make_unique<float []>(total_size);
//  ytot = std::make_unique<float []>(total_size);
//  ztot = std::make_unique<float []>(total_size);
//
//  bgtot = std::make_unique<float []>(total_size);
//  lltot = std::make_unique<float []>(total_size);
//  photot = std::make_unique<float []>(total_size);
//  crlbxtot = std::make_unique<float []>(total_size);
//  crlbytot = std::make_unique<float []>(total_size);

  printf("size:\t %d --- %d \n", total_size - xy_tmp_size, total_size);

  uint32_t j = total_size - xy_tmp_size;
  for (uint32_t i = 0; i < xy_tmp_size; i++, j++) {
    xtot[j] = xtmp[i];
    ytot[j] = ytmp[i];
    ztot[j] = z[i];
    //ztot[j] = z[i] + total_frames*num_dataset + 1;
    //ztot[j] = z[i] + FRAMES*num_dataset + 1;

    bgtot[j] = sub_background[i];
    lltot[j] = -2 * LL2[i];
    photot[j] = sub_photon[i];
    crlbxtot[j] = sqrt(CRLB2[i + 0 * new_size]);
    crlbytot[j] = sqrt(CRLB2[i + 1 * new_size]);
  }

  //printf("total_size: %d\n", total_size);

  //clean up
  //printf("clearing temp memory\n");

  free(xtmp);
  free(ytmp);
}

void SRscmos::saveData(const char* output_file_path) {
  ////// IMAGE GENERATION //////////
  uint32_t sz = this->input_data_width;
  uint32_t segnum = 64;

  double* rch = (double*)calloc(zm_color*sz * zm_color*sz, sizeof(double));
  double* gch = (double*)calloc(zm_color*sz * zm_color*sz, sizeof(double));
  double* bch = (double*)calloc(zm_color*sz * zm_color*sz, sizeof(double));

  printf("starting reconstruction\n");
  srhist_color(sz, zm_color, total_size, xtot, ytot, ztot, segnum, rch, gch, bch);
  printf("finished reconstruction\n");

  unsigned long im_size = sz * zm_color;

  dip::Image r(dip::NonOwnedRefToDataSegment((void *)rch), (void *)rch, dip::DT_DFLOAT, dip::UnsignedArray{ im_size, im_size });
  dip::Image g(dip::NonOwnedRefToDataSegment((void *)gch), (void *)gch, dip::DT_DFLOAT, dip::UnsignedArray{ im_size, im_size });
  dip::Image b(dip::NonOwnedRefToDataSegment((void *)bch), (void *)bch, dip::DT_DFLOAT, dip::UnsignedArray{ im_size, im_size });

  dip::Image rd, gd, bd;

  dip::Gauss(r, rd);
  dip::Gauss(g, gd);
  dip::Gauss(b, bd);

  double *rchsm = (double *)(rd.Data());
  double *gchsm = (double *)(gd.Data());
  double *bchsm = (double *)(bd.Data());

  double* rchsmst = imstretch_linear(im_size, rchsm, 0, 3, 0, 255);
  double* gchsmst = imstretch_linear(im_size, gchsm, 0, 3, 0, 255);
  double* bchsmst = imstretch_linear(im_size, bchsm, 0, 3, 0, 255);

  writeImage(output_file_path, rchsmst, gchsmst, bchsmst, im_size, im_size);
  //printf("size: %lu\n", im_size);

  free(rchsmst);
  free(gchsmst);
  free(bchsmst);

  free(rch);
  free(gch);
  free(bch);

  // useful for plotting
  //saveArray( "output/bgtot.dat", bgtot, total_size );
  //saveArray( "output/photot.dat", photot, total_size );
  //saveArray( "output/crlbxtot.dat", crlbxtot, total_size );
  //saveArray( "output/crlbytot.dat", crlbytot, total_size );
  //saveArray( "output/lltot.dat", lltot, total_size );

  //save2D( "output/xytot.dat", xtot, ytot, total_size );

}

void SRscmos::mapData(const SRscmos &local_data, uint32_t num_dataset, uint32_t frames_per_worker_fixed, uint32_t remaining_frames) {

  // find the lowres max 
  for (uint32_t i = 0; i < this->imsz; i++)
    for (uint32_t j = 0; j < this->imsz; j++) {
      if (this->maxprojim[i * this->imsz + j] < local_data.maxprojim[i * local_data.imsz + j]) {
        this->maxprojim[i * this->imsz + j] = local_data.maxprojim[i * local_data.imsz + j];
    }
  }

  for (uint32_t i = 0; i < this->imsz; i++) {
    for (uint32_t j = 0; j < this->imsz; j++) {
      this->sumim[i * this->imsz + j] = this->sumim[i * this->imsz + j] + local_data.sumim[i*local_data.imsz + j];
    }
  }

  for (uint32_t i = 0; i < this->imszzm; i++) {
    for (uint32_t j = 0; j < this->imszzm; j++) {
      this->imtot[i*this->imszzm + j] = this->imtot[i*this->imszzm + j] + local_data.imtot[i*local_data.imszzm + j];
    }
  }

  // append results to previous data
  this->total_size = this->total_size + local_data.total_size;

//  xtot = std::make_unique<float []>(total_size);
//  ytot = std::make_unique<float []>(total_size);
//  ztot = std::make_unique<float []>(total_size);
//
//  bgtot = std::make_unique<float []>(total_size);
//  lltot = std::make_unique<float []>(total_size);
//  photot = std::make_unique<float []>(total_size);
//  crlbxtot = std::make_unique<float []>(total_size);
//  crlbytot = std::make_unique<float []>(total_size);

//  float *xtot_ptr = this->xtot.get();
//  float *ytot_ptr = this->ytot.get();
//  float *ztot_ptr = this->ztot.get();
//
//  float *bgtot_ptr = this->bgtot.get();
//  float *lltot_ptr = this->lltot.get();
//  float *photot_ptr =   this->photot.get();
//  float *crlbxtot_ptr = this->crlbxtot.get();
//  float *crlbytot_ptr = this->crlbytot.get();
//
//  xtot_ptr = (float *)realloc(xtot_ptr, total_size * sizeof(float));
//  ytot_ptr = (float *)realloc(ytot_ptr, total_size * sizeof(float));
//  ztot_ptr = (float *)realloc(ztot_ptr, total_size * sizeof(float));
//
//  bgtot_ptr = (float *)realloc(bgtot_ptr, total_size * sizeof(float));
//  lltot_ptr = (float *)realloc(lltot_ptr, total_size * sizeof(float));
//  photot_ptr = (float *)realloc(photot_ptr, total_size * sizeof(float));
//  crlbxtot_ptr = (float *)realloc(crlbxtot_ptr, total_size * sizeof(float));
//  crlbytot_ptr = (float *)realloc(crlbytot_ptr, total_size * sizeof(float));

  xtot = (float *)realloc(xtot, total_size * sizeof(float));
  ytot = (float *)realloc(ytot, total_size * sizeof(float));
  ztot = (float *)realloc(ztot, total_size * sizeof(float));

  bgtot = (float *)realloc(bgtot, total_size * sizeof(float));
  lltot = (float *)realloc(lltot, total_size * sizeof(float));
  photot = (float *)realloc(photot, total_size * sizeof(float));
  crlbxtot = (float *)realloc(crlbxtot, total_size * sizeof(float));
  crlbytot = (float *)realloc(crlbytot, total_size * sizeof(float));

  printf("size:\t %d --- %d \n", total_size - local_data.total_size, total_size);

  int j = total_size - local_data.total_size;
  int correction = 1;

  if (num_dataset >= remaining_frames) {
    correction += remaining_frames;
  }

  for (uint32_t i = 0; i < local_data.total_size; i++, j++) {
    xtot[j] = local_data.xtot[i];
    ytot[j] = local_data.ytot[i];
    ztot[j] = (local_data.ztot[i] + frames_per_worker_fixed*num_dataset) + correction;

    bgtot[j] = local_data.bgtot[i];
    lltot[j] = local_data.lltot[i];
    photot[j] = local_data.photot[i];
    crlbxtot[j] = local_data.crlbxtot[i];
    crlbytot[j] = local_data.crlbytot[i];
  }

  printf("total_size: %d\n", total_size);
}

//void SRscmos::mergeData(float* cp_imtot, float* cp_xtot, float* cp_ytot, float* cp_ztot, float* cp_bgtot, float* cp_lltot, float* cp_photot, float* cp_crlbxtot, float* cp_crlbytot, double* cp_maxprojim, double* cp_sumim, uint32_t cp_total_size, uint32_t num_dataset, uint32_t total_frames, uint32_t remaining_frames, uint32_t data_width) {
//
//  /////
//  imsz = data_width - 3;
//  zm = 1 / (double)SRImage_pixelsize;
//  imszzm = (uint32_t)(imsz*zm);
//
//  imtot = std::make_unique<float []>(imszzm * imszzm);
//
//  maxprojim = std::make_unique<double[]>(imsz*imsz);
//  sumim = std::make_unique<double[]>(imsz*imsz);
//  /////
//
//  // find the lowres max 
//  for (uint32_t i = 0; i < imsz; i++)
//    for (uint32_t j = 0; j < imsz; j++) {
//      if (maxprojim[i * imsz + j] < cp_maxprojim[i * imsz + j]) {
//        maxprojim[i * imsz + j] = cp_maxprojim[i * imsz + j];
//    }
//  }
//
//  for (uint32_t i = 0; i < imsz; i++) {
//    for (uint32_t j = 0; j < imsz; j++) {
//      sumim[i * imsz + j] = sumim[i * imsz + j] + cp_sumim[i * imsz + j];
//    }
//  }
//
//  for (uint32_t i = 0; i < imszzm; i++) {
//    for (uint32_t j = 0; j < imszzm; j++) {
//      imtot[i*imszzm + j] = imtot[i*imszzm + j] + cp_imtot[i*imszzm + j];
//    }
//  }
//
//  // append results to previous data
//  total_size = total_size + cp_total_size;
//
////  xtot = std::make_unique<float []>(total_size);
////  ytot = std::make_unique<float []>(total_size);
////  ztot = std::make_unique<float []>(total_size);
////
////  bgtot = std::make_unique<float []>(total_size);
////  lltot = std::make_unique<float []>(total_size);
////  photot = std::make_unique<float []>(total_size);
////  crlbxtot = std::make_unique<float []>(total_size);
////  crlbytot = std::make_unique<float []>(total_size);
//
//  xtot = (float *)realloc(xtot, total_size * sizeof(float));
//  ytot = (float *)realloc(ytot, total_size * sizeof(float));
//  ztot = (float *)realloc(ztot, total_size * sizeof(float));
//
//  bgtot = (float *)realloc(bgtot, total_size * sizeof(float));
//  lltot = (float *)realloc(lltot, total_size * sizeof(float));
//  photot = (float *)realloc(photot, total_size * sizeof(float));
//  crlbxtot = (float *)realloc(crlbxtot, total_size * sizeof(float));
//  crlbytot = (float *)realloc(crlbytot, total_size * sizeof(float));
//
//  printf("size:\t %d --- %d \n", total_size - cp_total_size, total_size);
//
//  int j = total_size - cp_total_size;
//
//  int correction = 1;
//
//  if (num_dataset >= remaining_frames) {
//    correction += remaining_frames;
//  }
//
//  for (uint32_t i = 0; i < cp_total_size; i++, j++) {
//    xtot[j] = cp_xtot[i];
//    ytot[j] = cp_ytot[i];
//
//    ztot[j] = (cp_ztot[i] + total_frames*num_dataset) + correction;
//    //ztot[j] = cp_ztot[i] + total_frames*num_dataset + 1;
//    //ztot[j] = cp_ztot[i] + FRAMES*num_dataset + 1;
//
//    bgtot[j] = cp_bgtot[i];
//    lltot[j] = cp_lltot[i];
//
//    photot[j] = cp_photot[i];
//
//    crlbxtot[j] = cp_crlbxtot[i];
//    crlbytot[j] = cp_crlbytot[i];
//  }
//
//  printf("total_size: %d\n", total_size);
//}

//#define COMPRESSION_NONE 0

void SRscmos::writeMatFile(const std::string &gain_path, const std::string &var_path, const std::string &input_path, const std::string &output_path, uint32_t total_frames) {

  std::string output_mat_path = output_path;
  printf("saving to: %s\n", output_mat_path.c_str());

  mat_t *matfp;

  matvar_t *matvar, *field;
  size_t dims[2] = {(size_t)total_size,1};
  size_t struct_dims[2] = { 1,1 };

  const char *fieldnames[33] = { "datafolder", "resultfolder",
    "var_cal_file", "gain_cal_file",
    "xstart", "ystart",
    "imagesz", "startfrm", "totframes", "min_photon", "FitBoxSize",
    "iterations", "FitType", "PSFSigma", "SRImage_pixelsize",
    "llthreshold", "subregion_cropx", "subregion_cropy",
    "N_photon_min", "N_photon_max", "loc_uncer_max", "loc_uncer_min",
    "loc_x", "loc_y", "loc_z", "loc_bg", "loc_photons", "loc_uncerx",
    "loc_uncery", "loc_ll", "srim", "maxprojim", "sumim" };
  unsigned nfields = 33;

  matfp = Mat_CreateVer(output_mat_path.c_str(), NULL, MAT_FT_MAT73);

  if (NULL == matfp) {
    fprintf(stderr, "Error creating MAT file \"data.mat\"\n");
    //return EXIT_FAILURE;
    return;
  }

  matvar = Mat_VarCreateStruct("srobj", 2, struct_dims, fieldnames, nfields);
  if (NULL == matvar) {
    fprintf(stderr, "Error creating variable for a\n");
    Mat_Close(matfp);
    //return EXIT_FAILURE;
    return;
  }

  //datafolder
  size_t dim[2] = { 1, strlen(input_path.c_str()) };
  field = Mat_VarCreate("datafolder", MAT_C_CHAR, MAT_T_UTF8, 2, dim, (char*)(input_path.c_str()), 0);
  Mat_VarSetStructFieldByName(matvar, "datafolder", 0, field);

  //resultfolder
  dim[1] = strlen(output_path.c_str());
  field = Mat_VarCreate("resultfolder", MAT_C_CHAR, MAT_T_UTF8, 2, dim, (char*)(output_path.c_str()), 0);
  Mat_VarSetStructFieldByName(matvar, "resultfolder", 0, field);

  //var_cal_file
  dim[1] = strlen(var_path.c_str());
  field = Mat_VarCreate("var_cal_file", MAT_C_CHAR, MAT_T_UTF8, 2, dim, (char*)(var_path.c_str()), 0);
  Mat_VarSetStructFieldByName(matvar, "var_cal_file", 0, field);

  //gain_cal_file
  dim[1] = strlen(gain_path.c_str());
  field = Mat_VarCreate("gain_cal_file", MAT_C_CHAR, MAT_T_UTF8, 2, dim, (char*)(gain_path.c_str()), 0);
  Mat_VarSetStructFieldByName(matvar, "gain_cal_file", 0, field);

  // loc_x 
  field = Mat_VarCreate("loc_x", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, xtot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_x", 0, field);

  // loc_y
  field = Mat_VarCreate("loc_y", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, ytot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_y", 0, field);

  // loc_z 
  field = Mat_VarCreate("loc_z", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, ztot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_z", 0, field);

  // loc_bg 
  field = Mat_VarCreate("loc_bg", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, bgtot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_bg", 0, field);

  // loc_photons 
  field = Mat_VarCreate("loc_photons", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, photot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_photons", 0, field);

  // loc_uncerx 
  field = Mat_VarCreate("loc_uncerx", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, crlbxtot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_uncerx", 0, field);

  // loc_uncerx 
  field = Mat_VarCreate("loc_uncery", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, crlbytot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_uncery", 0, field);

  // loc_ll
  field = Mat_VarCreate("loc_ll", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, lltot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_ll", 0, field);

  //imagesz
  uint32_t imagesz = this->input_data_width;
  size_t dims2[2] = { 1,1 };
  field = Mat_VarCreate("imagesz", MAT_C_INT32, MAT_T_INT32, 2, dims2, &imagesz, 0);
  Mat_VarSetStructFieldByName(matvar, "imagesz", 0, field);

  // xstart
  uint32_t xstart = 1;
  field = Mat_VarCreate("xstart", MAT_C_INT32, MAT_T_INT32, 2, dims2, &xstart, 0);
  Mat_VarSetStructFieldByName(matvar, "xstart", 0, field);

  // ystart
  uint32_t ystart = 1;
  field = Mat_VarCreate("ystart", MAT_C_INT32, MAT_T_INT32, 2, dims2, &ystart, 0);
  Mat_VarSetStructFieldByName(matvar, "ystart", 0, field);

  // startfrm 
  uint32_t startfrm = 1;
  field = Mat_VarCreate("startfrm", MAT_C_INT32, MAT_T_INT32, 2, dims2, &startfrm, 0);
  Mat_VarSetStructFieldByName(matvar, "startfrm", 0, field);

  // totframes
  uint32_t totframes = total_frames;
  field = Mat_VarCreate("totframes", MAT_C_INT32, MAT_T_INT32, 2, dims2, &totframes, 0);
  Mat_VarSetStructFieldByName(matvar, "totframes", 0, field);

  // min_photon
  uint32_t min_photon_tmp = this->min_photon;
  field = Mat_VarCreate("min_photon", MAT_C_INT32, MAT_T_INT32, 2, dims2, &min_photon_tmp, 0);
  Mat_VarSetStructFieldByName(matvar, "min_photon", 0, field);

  // FitBoxSize
  uint32_t FitBoxSize = FIT_BOX_SIZE;
  field = Mat_VarCreate("FitBoxSize", MAT_C_INT32, MAT_T_INT32, 2, dims2, &FitBoxSize, 0);
  Mat_VarSetStructFieldByName(matvar, "FitBoxSize", 0, field);

  // iterations
  uint32_t iterations = ITERATIONS;
  field = Mat_VarCreate("iterations", MAT_C_INT32, MAT_T_INT32, 2, dims2, &iterations, 0);
  Mat_VarSetStructFieldByName(matvar, "iterations", 0, field);

  // FitType 
  uint32_t FitType = FIT_TYPE;
  field = Mat_VarCreate("FitType", MAT_C_INT32, MAT_T_INT32, 2, dims2, &FitType, 0);
  Mat_VarSetStructFieldByName(matvar, "FitType", 0, field);

  // PSFSigma
  field = Mat_VarCreate("PSFSigma", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims2, &(this->PSFSigma), 0);
  Mat_VarSetStructFieldByName(matvar, "PSFSigma", 0, field);

  // SRImage_pixelsize
  //float SRmage_pixelsize = 0.1000;
  field = Mat_VarCreate("SRImage_pixelsize", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims2, &(this->SRImage_pixelsize), 0);
  Mat_VarSetStructFieldByName(matvar, "SRImage_pixelsize", 0, field);

  // llthreshold
  field = Mat_VarCreate("llthreshold", MAT_C_INT32, MAT_T_INT32, 2, dims2, &(this->ll_threshold), 0);
  Mat_VarSetStructFieldByName(matvar, "llthreshold", 0, field);

  // subregion_cropx
  field = Mat_VarCreate("subregion_cropx", MAT_C_INT32, MAT_T_INT32, 2, dims2, &(this->subregion_cropx), 0);
  Mat_VarSetStructFieldByName(matvar, "subregion_cropx", 0, field);

  // subregion_cropy
  field = Mat_VarCreate("subregion_cropy", MAT_C_INT32, MAT_T_INT32, 2, dims2, &(this->subregion_cropy), 0);
  Mat_VarSetStructFieldByName(matvar, "subregion_cropy", 0, field);

  // N_photon_min
  field = Mat_VarCreate("N_photon_min", MAT_C_INT32, MAT_T_INT32, 2, dims2, &(this->N_photon_min), 0);
  Mat_VarSetStructFieldByName(matvar, "N_photon_min", 0, field);

  // N_photon_max
  field = Mat_VarCreate("N_photon_max", MAT_C_INT32, MAT_T_INT32, 2, dims2, &(this->N_photon_max), 0);
  Mat_VarSetStructFieldByName(matvar, "N_photon_max", 0, field);

  // loc_uncer_max 
  field = Mat_VarCreate("loc_uncer_max", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims2, &(this->loc_uncer_max), 0);
  Mat_VarSetStructFieldByName(matvar, "loc_uncer_max", 0, field);

  // loc_uncer_min 
  field = Mat_VarCreate("loc_uncer_min", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims2, &(this->loc_uncer_min), 0);
  Mat_VarSetStructFieldByName(matvar, "loc_uncer_min", 0, field);

  //srim
  //  float srim[imszzm][imszz];
  //  for( int i = 0; i < imszzm; i++ ) {
  //    for( int j = 0; j < imszzm; j++ ) {
  //      imtot[imszzm][imszzm] = imtot[i*imszzm+j];
  //    }
  //  }

  //srim
  //size_t dim2d[2] = {imszzm, imszzm};
  size_t dim2d[2];
  dim2d[0] = imszzm; dim2d[1] = imszzm;

  field = Mat_VarCreate("srim", MAT_C_SINGLE, MAT_T_SINGLE, 2, dim2d, imtot, 0);
  Mat_VarSetStructFieldByName(matvar, "srim", 0, field);

  //maxprojim
  dim2d[0] = this->input_data_width-3; dim2d[1] = this->input_data_width-3;
  field = Mat_VarCreate("maxprojim", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d, maxprojim, 0);
  Mat_VarSetStructFieldByName(matvar, "maxprojim", 0, field);

  //sumim
  //size_t dim2d[2] = {this->input_data_width,this->input_data_width};
  dim2d[0] = this->input_data_width-3; dim2d[1] = this->input_data_width-3;
  field = Mat_VarCreate("sumim", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d, sumim, 0);
  Mat_VarSetStructFieldByName(matvar, "sumim", 0, field);

  Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);
  Mat_VarFree(matvar);
  Mat_Close(matfp);
}


