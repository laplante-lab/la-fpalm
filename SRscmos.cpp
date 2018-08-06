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

  // set default parameters
  imsz = IMAGE_SIZE - 3; // 256 - 3 = 253

  zm = 1 / (double)SRImage_pixelsize;
  imszzm = (int)(imsz*zm);

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

  total_size = 0;
}

void SRscmos::setFitParameters(int min_photon, float PSFSigma, int zm_color, float SRImage_pixelsize) {
  this->min_photon = min_photon;
  this->PSFSigma = PSFSigma;
  this->zm_color = zm_color;
  this->SRImage_pixelsize = SRImage_pixelsize;

  zm = 1 / (double)SRImage_pixelsize;
  imszzm = (int)(imsz*zm);

  imtot = (float*)realloc( imtot, imszzm*imszzm*sizeof(float) );
}

void SRscmos::setThresholds( int ll_threshold, int subregion_cropx, int subregion_cropy, int N_photon_min, int N_photon_max, float loc_uncer_min, float loc_uncer_max ) {
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

// Use for reference
//void SRscmos::getCenterSubMatrix(Image<double> *input, Image<double> *output) {
//  uint32_t test_height = 16;
//  uint32_t test_width = 16;
//
//  Image<uint32_t> *test = new Image<uint32_t>(1,test_height,test_width);
//  for(int i = 0; i < test_height; i++) {
//    for(int j = 0; j < test_width; j++) {
//      (*test)[{i,j}] = (i*test_width)+j;
//    }
//  }
//  std::cout << (*test) << std::endl;
//
//  // 4x4 center array
//
//  int is_corner = test_width%2==0;
//  int center_index_x = ceil(test_width/2) - is_corner;
//  int center_index_y = ceil(test_width/2) - is_corner;
//
//  printf("%d, center: [%d, %d]\n", is_corner, center_index_x, center_index_y);
//  printf("center: %d\n", (*test)[{center_index_x, center_index_y}]); 
//  
//  uint32_t center_height = 2;
//  uint32_t center_width = 2;
//  Image<uint32_t> *center = new Image<uint32_t>(1, center_height, center_width); 
//
//  if ( center_height == 1 && center_width == 1 ) {
//      (*center)[{0,0}] = (*test)[{center_index_x, center_index_y}]; 
//  }
//  else {
//    for( uint32_t i = 0; i < center->height; i++ ) {
//      for( uint32_t j = 0; j < center->height; j++ ) {
//        (*center)[{i,j}] = (*test)[{i+ceil((double)test_height/center_height),j+ceil((double)test_width/center_width)}];
//      }
//    }
//  }
//
//  std::cout << (*center) << std::endl;
//
//  for ( uint32_t i = 0; i < 5; i++ ) {
//    //printf("%lf\n", (*fin)[std::make_tuple(0,i,0)]);
//    printf("%lf %lf\n", (*nm->o_map)[{i+128,0+128}], (*nm->gain_map)[{i+128,0+128}]);
//    //printf("%lf %lf\n", (*nm->o_map)[{i+mapbase,0+mapbase}], (*nm->gain_map)[{i+mapbase,0+mapbase}]);
//  }
//}

void SRscmos::run(NoiseMap* nm, Image<uint16_t> *frame_stack, int num_dataset) {
  // This applies the noise filter to the input images
  Image<double> *fin = new Image<double>( frame_stack->num_frames, frame_stack->height, frame_stack->width ); 

  int mapbase = 0;

  // crop out subsection of camera noise
  if ( frame_stack->height < nm->gain_map->height ) {
    mapbase = (nm->gain_map->height - frame_stack->height)/2;
  }

  printf("mapbase: %d\n", mapbase);

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

  printf("convolution done\n");

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

  printf("dilation done\n");

  //////////////////
  double *a = (double *)malloc(sizeof(double));
  
  int num_regions = 0;

  for (uint32_t it = 0; it < im_max->num_frames; it++)
    for (uint32_t j = 0; j < im_max->height; j++)
      for (uint32_t i = 0; i < im_max->width; i++)
        if ((*im_max)[std::make_tuple(it,i,j)]) {
          num_regions++;
          a = (double *)realloc(a, num_regions * sizeof(double));
          a[num_regions - 1] = to1D(i, j, it, im_max->height, im_max->width);
        }

  double *x = (double *)malloc(num_regions * sizeof(double));
  double *y = (double *)malloc(num_regions * sizeof(double));
  double *z = (double *)malloc(num_regions * sizeof(double));

  uint32_t imsz = im_max->width;
  delete im_max;

  printf("num_regions: %d\n", num_regions);

  for (int i = 0; i < num_regions; i++) {
    z[i] = floor(a[i] / 253 / 253);
    x[i] = fmod(fmod(a[i], (imsz*imsz)), imsz);
    y[i] = floor(fmod(a[i], (imsz*imsz)) / imsz);
    //printf("%0.4lf  %0.4lf  %0.4lf\n", x[i], y[i], z[i]);
  }

  free(a);

  //////////////////
  // sidemask //
  bool *sidemask = (bool *)malloc(num_regions * sizeof(bool));
  int new_num_regions = num_regions;

  for (int i = 0; i < num_regions; i++) {
    sidemask[i] = (x[i] == 0) || (x[i] == frame_stack->width - 2 - 1) || (y[i] == 0) || (y[i] == frame_stack->height - 2 - 1);
    if (sidemask[i]) {
      new_num_regions--;
    }
  }

  double *temp_x = (double *)malloc(new_num_regions * sizeof(double));
  double *temp_y = (double *)malloc(new_num_regions * sizeof(double));
  double *temp_z = (double *)malloc(new_num_regions * sizeof(double));

  int ind = 0;
  for (int i = 0; i < num_regions; i++) {
    if (!sidemask[i]) {
      temp_x[ind] = x[i];
      temp_y[ind] = y[i];
      temp_z[ind++] = z[i];
      //printf("%3lf %3lf %3lf\n", x[i], y[i], z[i]);
    }
  }

  free(sidemask);

  free(x);
  free(y);
  free(z);

  x = temp_x;
  y = temp_y;
  z = temp_z;
  num_regions = new_num_regions;
  ////////////////////////////

  //printf("num_regions: %d\n", num_regions);

  ///////// isolate sub-regions /////////////
  //int boxsz = 7;

  printf("first: x[0]:%0.1lf    y[0]:%0.1lf    z[0]:%0.1lf\n", x[0], y[0], z[0] );
  printf("center: %0.4lf\n", (*fin)[std::make_tuple((int)z[0],(int)x[0],(int)y[0])] );
  printf("num_regions:%d\n", num_regions);
  //exit(0);
  

  float* t2 = (float *)malloc(num_regions * sizeof(float));
  float* l2 = (float *)malloc(num_regions * sizeof(float));
  
  Image<float> *varmap_cut = new Image<float>(fcin->num_frames, fcin->height, fcin->width);
  Image<float> *gainmap_cut = new Image<float>(fcin->num_frames, fcin->height, fcin->width);

    //float *varmap_cut = (float*)malloc(253 * 253 * fcin->num_frames * sizeof(float));
  //  float *gainmap_cut = (float*)malloc(253 * 253 * total_frames * sizeof(float));
  // fcin = 253 x 253 x 1000, fin = 1000 x 253 x 253

  for (uint32_t i = 0; i < fcin->num_frames; i++)
    for (uint32_t j = 0; j < fcin->height; j++)
      for (uint32_t k = 0; k < fcin->width; k++) {
        (*varmap_cut)[std::make_tuple(i,k,j)] = (float)(*(nm->var_map))[{k+mapbase+1,j+mapbase+1}];
        (*gainmap_cut)[std::make_tuple(i,k,j)] = (float)(*(nm->gain_map))[{k+mapbase+1,j+mapbase+1}];
      }


  printf("check2\n");
  float *subregion = cMakeSubregions(x, y, z, num_regions, FIT_BOX_SIZE, fcin->pixels, l2, t2, fcin->height, fcin->width, fcin->num_frames );
  printf("check3\n");

  float* tmp_l = (float *)malloc(num_regions * sizeof(float));
  float* tmp_t = (float *)malloc(num_regions * sizeof(float));
 
  float *subvarim = cMakeSubregions(x, y, z, num_regions, FIT_BOX_SIZE, varmap_cut->pixels, tmp_l, tmp_t, fcin->height, fcin->width, fcin->num_frames);
  float *subgainim = cMakeSubregions(x, y, z, num_regions, FIT_BOX_SIZE, gainmap_cut->pixels, tmp_l, tmp_t, fcin->height, fcin->width, fcin->num_frames);

  delete varmap_cut;
  delete gainmap_cut;

  free(tmp_l); free(tmp_t);

  free(x);
  free(y);
  printf("checkall\n");

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

  float *sub_x = (float *)malloc(num_regions * sizeof(float));
  float *sub_y = (float *)malloc(num_regions * sizeof(float));
  float *sub_photon = (float *)malloc(num_regions * sizeof(float));
  float *sub_background = (float *)malloc(num_regions * sizeof(float));
  float *CRLB2 = (float *)malloc(4 * num_regions * sizeof(float));
  float *LL2 = (float *)malloc(num_regions * sizeof(float));

  printf("PSFSigma: %0.4lf, num_dataset: %d \n", PSFSigma, num_dataset);

  SRsCMOS_MLE(num_dataset, num_regions, (float *)subregion, (float *)subvarim, (float *)subgainim, FIT_TYPE, PSFSigma, ITERATIONS, sub_x, sub_y, sub_photon, sub_background, CRLB2, LL2);

  free(subregion);
  free(subvarim);
  free(subgainim);


  printf("finished GPU, %d\n", num_dataset);

  bool* llmask2 = (bool *)malloc(num_regions * sizeof(bool));
  bool* xymask = (bool *)malloc(num_regions * sizeof(bool));
  bool* intmask = (bool *)malloc(num_regions * sizeof(bool));

  for (int i = 0; i < num_regions; i++)
    llmask2[i] = (LL2[i] * (-2)) > ll_threshold;

  for (int i = 0; i < num_regions; i++)
    xymask[i] = (sub_x[i] >(FIT_BOX_SIZE - subregion_cropx - 1)) || (sub_x[i] < subregion_cropx) || (sub_y[i] > (FIT_BOX_SIZE - subregion_cropx - 1)) || (sub_y[i] < subregion_cropx);

  for (int i = 0; i < num_regions; i++)
    intmask[i] = (sub_photon[i] < N_photon_min) || (sub_photon[i] > N_photon_max);

  bool* uncermask = (bool *)malloc(num_regions * sizeof(bool));

  for (int i = 0; i < num_regions; i++)
    uncermask[i] = (sqrt(CRLB2[i + 1 * num_regions]) < loc_uncer_min) || (sqrt(CRLB2[i + 0 * num_regions]) > loc_uncer_max);

  bool* totmask = (bool *)malloc(num_regions * sizeof(bool));

  int mask_size = 0;
  for (int i = 0; i < num_regions; i++) {
    totmask[i] = llmask2[i] || xymask[i] || intmask[i] || uncermask[i];

    if (!totmask[i])
      mask_size++;
  }

  free(llmask2);
  free(xymask);
  free(intmask);
  free(uncermask);

  int new_size = mask_size;
  printf("new_size: %d\n", new_size);

  float* mask_sub_x = (float *)malloc(new_size * sizeof(float));
  float* mask_sub_y = (float *)malloc(new_size * sizeof(float));
  float* mask_sub_photon = (float *)malloc(new_size * sizeof(float));
  float* mask_sub_background = (float *)malloc(new_size * sizeof(float));
  float* mask_CRLB2 = (float *)malloc(4 * new_size * sizeof(float));
  float* mask_l2 = (float *)malloc(new_size * sizeof(float));
  float* mask_t2 = (float *)malloc(new_size * sizeof(float));
  double* mask_z = (double *)malloc(new_size * sizeof(double));
  float* mask_LL2 = (float *)malloc(new_size * sizeof(float));

  int indexer = 0;
  for (int i = 0; i < num_regions; i++)
    if (!totmask[i]) {
      mask_sub_x[indexer] = sub_x[i];
      mask_sub_y[indexer] = sub_y[i];
      mask_sub_photon[indexer] = sub_photon[i];
      mask_sub_background[indexer] = sub_background[i];
      for (int j = 0; j < 4; j++)
        mask_CRLB2[indexer + j*new_size] = CRLB2[i + j*num_regions];
      mask_l2[indexer] = (float)l2[i];
      mask_t2[indexer] = t2[i];
      mask_z[indexer] = z[i];
      mask_LL2[indexer] = LL2[i];

      indexer++;
    }
    
  free(totmask);
  //printf("checkpoint: %d\n", indexer);
  free(sub_x);
  free(sub_y);
  free(sub_photon);
  free(sub_background);
  free(CRLB2);
  free(l2);
  free(t2);
  free(z);
  free(LL2);

  sub_x = mask_sub_x;
  sub_y = mask_sub_y;
  sub_photon = mask_sub_photon;
  sub_background = mask_sub_background;
  CRLB2 = mask_CRLB2;
  l2 = mask_l2;
  t2 = mask_t2;
  z = mask_z;
  LL2 = mask_LL2;
  
  //printf("start reconstruction\n");
  //// Reconstruction /////
  double maxblob = 100000;
  int maxk = (int)ceil(new_size / maxblob);
  //printf("maxk:%d\n", maxk);
 
  ////////CONTINUE/////////
  float *xresult2, *yresult2;

  float* xtmp = (float*)malloc(sizeof(float));
  float* ytmp = (float*)malloc(sizeof(float));
  int xy_tmp_size = 0;

  /////// Parse data /////////////////
  for (int ii = 0; ii < maxk; ii++) {
    int bst = (ii)*maxblob + 1;
    //printf("bst: %d\n", bst);

    int bed;

    if (ii == (maxk - 1))
      bed = new_size;
    else
      bed = ii * maxblob;

    xresult2 = (float *)malloc((bed - bst + 1) * sizeof(float));
    yresult2 = (float *)malloc((bed - bst + 1) * sizeof(float));
 
    int indx = 0;
    for (int i = bst - 1; i < bed; i++) {
      xresult2[indx] = sub_y[i] + l2[i];
      yresult2[indx++] = sub_x[i] + t2[i];
    }

    //printf("bed-bst: %d\n", (bed-bst+1));
    //printf("indx: %d\n", indx);
    //printf("imszzm: %d\n", imszzm);

    float* tempx = (float*)malloc((bed - bst + 1) * sizeof(float));
    float* tempy = (float*)malloc((bed - bst + 1) * sizeof(float));
  
    for (int i = 0; i < indx; i++) {
      tempx[i] = xresult2[i] * zm;
      tempy[i] = yresult2[i] * zm;
    }
  
    xy_tmp_size = xy_tmp_size + indx;
    xtmp = (float*)realloc(xtmp, xy_tmp_size * sizeof(float));
    ytmp = (float*)realloc(ytmp, xy_tmp_size * sizeof(float));
  
    uint16_t *im2 = cHistRecon(imszzm, imszzm, indx, tempx, tempy, 1);
 
    for (int i = 0; i < imszzm; i++)
      for (int j = 0; j < imszzm; j++)
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

  xtot = (float *)realloc(xtot, total_size * sizeof(float));
  ytot = (float *)realloc(ytot, total_size * sizeof(float));
  ztot = (float *)realloc(ztot, total_size * sizeof(float));

  bgtot = (float *)realloc(bgtot, total_size * sizeof(float));
  lltot = (float *)realloc(lltot, total_size * sizeof(float));
  photot = (float *)realloc(photot, total_size * sizeof(float));
  crlbxtot = (float *)realloc(crlbxtot, total_size * sizeof(float));
  crlbytot = (float *)realloc(crlbytot, total_size * sizeof(float));

  //printf("size:\t %d --- %d \n", total_size - xy_tmp_size, total_size);

  int j = total_size - xy_tmp_size;
  for (int i = 0; i < xy_tmp_size; i++, j++) {
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
  //exit(0);

  //clean up
  //printf("clearing temp memory\n");

  free(xtmp);
  free(ytmp);

  free(sub_x);
  free(sub_y);
  free(sub_photon);
  free(sub_background);
  free(CRLB2);
  free(l2);
  free(t2);
  free(z);
  free(LL2);

}

void SRscmos::saveData(const char* output_file_path) {
  ////// IMAGE GENERATION //////////
  int sz = this->input_data_width;
  int segnum = 64;

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

  // useful for plotting
  //saveArray( "output/bgtot.dat", bgtot, total_size );
  //saveArray( "output/photot.dat", photot, total_size );
  //saveArray( "output/crlbxtot.dat", crlbxtot, total_size );
  //saveArray( "output/crlbytot.dat", crlbytot, total_size );
  //saveArray( "output/lltot.dat", lltot, total_size );

  //save2D( "output/xytot.dat", xtot, ytot, total_size );

}

void SRscmos::mergeData(float* cp_imtot, float* cp_xtot, float* cp_ytot, float* cp_ztot, float* cp_bgtot, float* cp_lltot, float* cp_photot, float* cp_crlbxtot, float* cp_crlbytot, double* cp_maxprojim, double* cp_sumim, int cp_total_size, int num_dataset, int total_frames, int remaining_frames) {

  // find the lowres max 
  for (int i = 0; i < 253; i++)
    for (int j = 0; j < 253; j++) {
      if (maxprojim[i * 253 + j] < cp_maxprojim[i * 253 + j]) {
        maxprojim[i * 253 + j] = cp_maxprojim[i * 253 + j];
    }
  }

  for (int i = 0; i < 253; i++) {
    for (int j = 0; j < 253; j++) {
      sumim[i * 253 + j] = sumim[i * 253 + j] + cp_sumim[i * 253 + j];
    }
  }

  for (int i = 0; i < imszzm; i++) {
    for (int j = 0; j < imszzm; j++) {
      imtot[i*imszzm + j] = imtot[i*imszzm + j] + cp_imtot[i*imszzm + j];
    }
  }

  // append results to previous data
  total_size = total_size + cp_total_size;

  xtot = (float *)realloc(xtot, total_size * sizeof(float));
  ytot = (float *)realloc(ytot, total_size * sizeof(float));
  ztot = (float *)realloc(ztot, total_size * sizeof(float));

  bgtot = (float *)realloc(bgtot, total_size * sizeof(float));
  lltot = (float *)realloc(lltot, total_size * sizeof(float));
  photot = (float *)realloc(photot, total_size * sizeof(float));
  crlbxtot = (float *)realloc(crlbxtot, total_size * sizeof(float));
  crlbytot = (float *)realloc(crlbytot, total_size * sizeof(float));

  printf("size:\t %d --- %d \n", total_size - cp_total_size, total_size);

  int j = total_size - cp_total_size;

  int correction = 1;

  if (num_dataset >= remaining_frames) {
    correction += remaining_frames;
  }

  for (int i = 0; i < cp_total_size; i++, j++) {
    xtot[j] = cp_xtot[i];
    ytot[j] = cp_ytot[i];

    ztot[j] = (cp_ztot[i] + total_frames*num_dataset) + correction;
    //ztot[j] = cp_ztot[i] + total_frames*num_dataset + 1;
    //ztot[j] = cp_ztot[i] + FRAMES*num_dataset + 1;

    bgtot[j] = cp_bgtot[i];
    lltot[j] = cp_lltot[i];

    photot[j] = cp_photot[i];

    crlbxtot[j] = cp_crlbxtot[i];
    crlbytot[j] = cp_crlbytot[i];
  }

  printf("total_size: %d\n", total_size);
}

//#define COMPRESSION_NONE 0

void SRscmos::writeMatFile(const std::string &gain_path, const std::string &var_path, const std::string &input_path, const std::string &output_path, int total_frames) {

  std::string output_mat_path = output_path;
  printf("saving to: %s\n", output_mat_path.c_str());
//  output_mat_path += "data.mat";

  mat_t *matfp;

  matvar_t *matvar, *field;
  //size_t dims[2] = {total_size,1};
  size_t dims[2];
  dims[0] = total_size; dims[1] = 1;

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

  //matfp = Mat_CreateVer(output_mat_path.c_str(), NULL, MAT_FT_DEFAULT);
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
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, xtot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_x", 0, field);

  // loc_y
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, ytot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_y", 0, field);

  // loc_z 
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, ztot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_z", 0, field);

  // loc_bg 
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, bgtot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_bg", 0, field);

  // loc_photons 
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, photot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_photons", 0, field);

  // loc_uncerx 
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, crlbxtot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_uncerx", 0, field);

  // loc_uncerx 
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, crlbytot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_uncery", 0, field);

  // loc_ll
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, lltot, 0);
  Mat_VarSetStructFieldByName(matvar, "loc_ll", 0, field);

  //imagesz
  int imagesz = this->input_data_width;
  size_t dims2[2] = { 1,1 };
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &imagesz, 0);
  Mat_VarSetStructFieldByName(matvar, "imagesz", 0, field);

  // xstart
  int xstart = 1;
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &xstart, 0);
  Mat_VarSetStructFieldByName(matvar, "xstart", 0, field);

  // ystart
  int ystart = 1;
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &ystart, 0);
  Mat_VarSetStructFieldByName(matvar, "ystart", 0, field);

  // startfrm 
  int startfrm = 1;
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &startfrm, 0);
  Mat_VarSetStructFieldByName(matvar, "startfrm", 0, field);

  // totframes
  int totframes = total_frames;
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &totframes, 0);
  Mat_VarSetStructFieldByName(matvar, "totframes", 0, field);

  // min_photon
  int min_photon_tmp = this->min_photon;
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &min_photon_tmp, 0);
  Mat_VarSetStructFieldByName(matvar, "min_photon", 0, field);

  // FitBoxSize
  int FitBoxSize = FIT_BOX_SIZE;
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &FitBoxSize, 0);
  Mat_VarSetStructFieldByName(matvar, "FitBoxSize", 0, field);

  // iterations
  int iterations = ITERATIONS;
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &iterations, 0);
  Mat_VarSetStructFieldByName(matvar, "iterations", 0, field);

  // FitType 
  int FitType = FIT_TYPE;
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &FitType, 0);
  Mat_VarSetStructFieldByName(matvar, "FitType", 0, field);

  // PSFSigma
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims2, &(this->PSFSigma), 0);
  Mat_VarSetStructFieldByName(matvar, "PSFSigma", 0, field);

  // SRImage_pixelsize
  //float SRmage_pixelsize = 0.1000;
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims2, &(this->SRImage_pixelsize), 0);
  Mat_VarSetStructFieldByName(matvar, "SRImage_pixelsize", 0, field);

  // llthreshold
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &(this->ll_threshold), 0);
  Mat_VarSetStructFieldByName(matvar, "llthreshold", 0, field);

  // subregion_cropx
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &(this->subregion_cropx), 0);
  Mat_VarSetStructFieldByName(matvar, "subregion_cropx", 0, field);

  // subregion_cropy
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &(this->subregion_cropy), 0);
  Mat_VarSetStructFieldByName(matvar, "subregion_cropy", 0, field);

  // N_photon_min
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &(this->N_photon_min), 0);
  Mat_VarSetStructFieldByName(matvar, "N_photon_min", 0, field);

  // N_photon_max
  field = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims2, &(this->N_photon_max), 0);
  Mat_VarSetStructFieldByName(matvar, "N_photon_max", 0, field);

  // loc_uncer_max 
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims2, &(this->loc_uncer_max), 0);
  Mat_VarSetStructFieldByName(matvar, "loc_uncer_max", 0, field);

  // loc_uncer_min 
  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims2, &(this->loc_uncer_min), 0);
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

  field = Mat_VarCreate(NULL, MAT_C_SINGLE, MAT_T_SINGLE, 2, dim2d, imtot, 0);
  Mat_VarSetStructFieldByName(matvar, "srim", 0, field);

  //maxprojim
  dim2d[0] = 253; dim2d[1] = 253;
  field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d, maxprojim, 0);
  Mat_VarSetStructFieldByName(matvar, "maxprojim", 0, field);

  //sumim
  //size_t dim2d[2] = {253,253};
  dim2d[0] = 253; dim2d[1] = 253;
  field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d, sumim, 0);
  Mat_VarSetStructFieldByName(matvar, "sumim", 0, field);


  Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);
  //Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);
  Mat_VarFree(matvar);
  Mat_Close(matfp);
  return;
  //return EXIT_SUCCESS;
}


