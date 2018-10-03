/**
*  @file    ComputeManager.cpp
*  @author  John Ravi (jjravi)
*  @date    6/12/2018
*  @version 2.0 
*
*  @brief Manages multiple independent streams of FPALM analysis. 
*
*  @section DESCRIPTION
*
*  Allows for different modes of execution for the analysis. 
*/

#include "ComputeManager.h"

ComputeManager::ComputeManager() {
  input = NULL;
}

ComputeManager::~ComputeManager() {
}

void ComputeManager::compute(const std::string &data_path, const std::string &output_path) {

  const std::string &inversegain_path = "./noise_map/la-gainmap.mat";
  const std::string &varmap_path = "./noise_map/la-varmap.mat";

  /// Step 1: sCMOS camera noise ///
  NoiseMap *nm = new NoiseMap();
  nm->parseDataMat( inversegain_path, varmap_path );

  /// Step 2: Load Images ///
  if ( !input ) {
    input = new DataIO(data_path);
  }

  printf("Data contains: %d frames\n", input->number_of_frames);
  printf("Frame size: %dx%d\n", input->frame_height, input->frame_width);

  int number_of_iterations = std::ceil((float)input->number_of_frames/MAX_FRAMES_PER_ITERATION);
  number_of_iterations = std::max(MIN_THREAD_COUNT, number_of_iterations);

#ifdef OMP
  omp_set_num_threads(number_of_iterations);
  printf("Launching %d threads.\n", number_of_iterations);
#endif

  int frames_per_worker = floor( input->number_of_frames/number_of_iterations );
  int remaining_frames = input->number_of_frames % number_of_iterations;

  printf("Analyzing %d frames per worker.\n", frames_per_worker);
  printf("Remaining %d frames, will be added to the first %d workers.\n", remaining_frames, remaining_frames);

  SRscmos *thread_sr[number_of_iterations];

#ifdef OMP
#pragma omp parallel for default(shared)
#endif
  for ( int num_dataset = 0; num_dataset < number_of_iterations; num_dataset++ ) {
    printf("\n------worker: %d-------\n", num_dataset);

    thread_sr[num_dataset] = new SRscmos();

    int offset = 0;
    int frames_per_worker_fixed = frames_per_worker;
    if ( num_dataset < remaining_frames )
      frames_per_worker_fixed++;
    else
      offset = remaining_frames;

    int start = (frames_per_worker_fixed*num_dataset + offset); 
    int end = (frames_per_worker_fixed*num_dataset + offset + frames_per_worker_fixed ); 
    int sub_num_frames = end - start;

    Image<uint16_t> *frame_stack;

#ifdef OMP
#pragma omp critical // faster to read sequentially
#endif
    { 
      printf("(%2d) reading frame_range: %6d -- %6d [%d]\n",num_dataset, start, end, sub_num_frames);
      frame_stack = input->parseImageData( start, end );
    }

    printf("Analyzing %d images (%dx%d) from %s\n", frame_stack->num_frames, frame_stack->height, frame_stack->width, input->data_path.c_str()); 

    thread_sr[num_dataset]->run(nm, frame_stack, num_dataset);
    if(frame_stack) delete frame_stack;
  }

#ifdef OMP
#pragma omp barrier
#endif

  SRscmos *sr = new SRscmos(*(thread_sr[0]));

  // merge results
  for ( int i = 0; i < number_of_iterations; i++ ) {
    SRscmos *one = thread_sr[i];

    int frames_per_worker_fixed = frames_per_worker;
    if ( i < remaining_frames )
      frames_per_worker_fixed++;

    sr->mapData( *(thread_sr[i]), i, frames_per_worker_fixed, remaining_frames);
    delete one;
  }

  if( createOutput( output_path ) != 0 ) {
    exit(EXIT_FAILURE);
  }

  ///// SAVE MATLAB FILE /////
  std::string output_mat_path = output_path;
  output_mat_path += "/data.mat";

  sr->writeMatFile( inversegain_path, varmap_path, data_path, output_mat_path, input->number_of_frames );
  printf("saved data in %s file\n", output_mat_path.c_str());

  if ( sr->total_size == 0 ) {
    printf("nothing to plot\n");
    return;
  }

  std::string output_image_path = output_path;
  output_image_path += "out.tiff";

  sr->saveData( output_image_path.c_str() );
  printf("outputted reconstructed image \n");

  printf("plotting results (multithreaded)\n");
  plot( output_path.c_str(), sr); 

  delete nm;
  delete sr;

  delete input;
}

void ComputeManager::compute(const std::string inversegain_path, const std::string varmap_path, const std::string data_path, const std::string output_path, 
                             int min_photon, float PSFSigma, int zm_color, float SRImage_pixelsize, int ll_threshold, 
                             int subregion_cropx, int subregion_cropy, int N_photon_min, int N_photon_max, float loc_uncer_min, float loc_uncer_max) {

  /// Step 1: sCMOS camera noise /////
  NoiseMap *nm = new NoiseMap();
  nm->parseDataMat( inversegain_path, varmap_path );

  ///// Step 2: Load Images //////
  if ( !input )
    input = new DataIO(data_path); // will know if it is a file or folder, and the number of frames

  printf("Data contains: %d frames\n", input->number_of_frames);
  printf("Frame size: %dx%d\n", input->frame_height, input->frame_width);

  int number_of_iterations = std::ceil((float)input->number_of_frames/MAX_FRAMES_PER_ITERATION);
  number_of_iterations = std::max(MIN_THREAD_COUNT, number_of_iterations);

#ifdef OMP
  omp_set_num_threads(number_of_iterations);
  printf("Launching %d threads.\n", number_of_iterations);
#endif

  int frames_per_worker = floor( input->number_of_frames/number_of_iterations );
  int remaining_frames = input->number_of_frames % number_of_iterations;

  printf("Analyzing %d frames per worker.\n", frames_per_worker);
  printf("Remaining %d frames, will be added to the first %d workers.\n", remaining_frames, remaining_frames);

  SRscmos *thread_sr[number_of_iterations];

#ifdef OMP
#pragma omp parallel for default(shared)
#endif
  for ( int num_dataset = 0; num_dataset < number_of_iterations; num_dataset++ ) {
    printf("\n------worker: %d-------\n", num_dataset);

    thread_sr[num_dataset] = new SRscmos();
    thread_sr[num_dataset]->setFitParameters( min_photon, PSFSigma, zm_color, SRImage_pixelsize);
    thread_sr[num_dataset]->setThresholds( ll_threshold, subregion_cropx, subregion_cropy, N_photon_min, N_photon_max, loc_uncer_min, loc_uncer_max);

    int offset = 0;
    int frames_per_worker_fixed = frames_per_worker;
    if ( num_dataset < remaining_frames )
      frames_per_worker_fixed++;
    else
      offset = remaining_frames;

    int start = (frames_per_worker_fixed*num_dataset + offset); 
    int end = (frames_per_worker_fixed*num_dataset + offset + frames_per_worker_fixed ); 
    int total_num_frames = end - start;

    Image<uint16_t> *frame_stack;

#ifdef OMP
#pragma omp critical
#endif
    { 
      printf("(%2d) reading frame_range: %6d -- %6d [%d]\n",num_dataset, start, end, total_num_frames);
      frame_stack = input->parseImageData( start, end );
    }

    printf("Analyzing %d images (%dx%d) from %s\n", frame_stack->num_frames, frame_stack->height, frame_stack->width, input->data_path.c_str()); 

    thread_sr[num_dataset]->run(nm, frame_stack, num_dataset);
    if(frame_stack) delete frame_stack;
  }

#ifdef OMP
#pragma omp barrier
#endif

  SRscmos *sr = new SRscmos(*(thread_sr[0]));
  sr->setFitParameters( min_photon, PSFSigma, zm_color, SRImage_pixelsize);
  sr->setThresholds( ll_threshold, subregion_cropx, subregion_cropy, N_photon_min, N_photon_max, loc_uncer_min, loc_uncer_max);

  // merge results
  for ( int i = 0; i < number_of_iterations; i++ ) {
    SRscmos *one = thread_sr[i];

    int frames_per_worker_fixed = frames_per_worker;
    if ( i < remaining_frames )
      frames_per_worker_fixed++;

    sr->mapData( *(thread_sr[i]), i, frames_per_worker_fixed, remaining_frames);
    delete one;
  }

  if( createOutput( output_path ) != 0 ) {
    exit(EXIT_FAILURE);
  }

  ///// SAVE MATLAB FILE /////
  std::string output_mat_path = output_path;
  output_mat_path += "/data.mat";

  sr->writeMatFile( inversegain_path, varmap_path, data_path, output_mat_path, input->number_of_frames );
  printf("saved data in %s file\n", output_mat_path.c_str());

  if ( sr->total_size == 0 ) {
    printf("nothing to plot\n");
    return;
  }

  std::string output_image_path = output_path;
  output_image_path += "out.tiff";

  sr->saveData( output_image_path.c_str() );
  printf("outputted reconstructed image \n");

  printf("plotting results (multithreaded)\n");
  plot( output_path.c_str(), sr); 

  delete nm;
  delete sr;

  delete input;
}

void ComputeManager::compute2(const std::string inversegain_path, const std::string varmap_path, const std::string data_path, const std::string output_path, 
                              int min_photon, float PSFSigma, int zm_color, float SRImage_pixelsize, int ll_threshold, 
                              int subregion_cropx, int subregion_cropy, int N_photon_min, int N_photon_max, float loc_uncer_min, float loc_uncer_max, 
                              int new_start, int new_end) {

  /// Step 1: sCMOS camera noise /////
  NoiseMap *nm = new NoiseMap();
  nm->parseDataMat( inversegain_path, varmap_path );

  ///// Step 2: Load Images //////
  if ( !input )
    input = new DataIO(data_path); // will know if it is a file or folder, and the number of frames

  input->number_of_frames = new_end - new_start;
  printf("Data contains: %d frames\n", input->number_of_frames);
  printf("Frame size: %dx%d\n", input->frame_height, input->frame_width);

  int number_of_iterations = std::ceil((float)input->number_of_frames/MAX_FRAMES_PER_ITERATION);
  number_of_iterations = std::max(MIN_THREAD_COUNT, number_of_iterations);

  if ( input->number_of_frames < (uint32_t)number_of_iterations )
    number_of_iterations = input->number_of_frames;

#ifdef OMP
  omp_set_num_threads(number_of_iterations);
  printf("Launching %d threads.\n", number_of_iterations);
#endif

  int frames_per_worker = floor( input->number_of_frames/number_of_iterations );
  int remaining_frames = input->number_of_frames % number_of_iterations;

  printf("Analyzing %d frames per worker.\n", frames_per_worker);
  printf("Remaining %d frames, will be added to the first %d workers.\n", remaining_frames, remaining_frames);

  SRscmos *thread_sr[number_of_iterations];

#ifdef OMP
#pragma omp parallel for default(shared)
#endif
  for ( int num_dataset = 0; num_dataset < number_of_iterations; num_dataset++ ) {
    printf("\n------worker: %d-------\n", num_dataset);

    thread_sr[num_dataset] = new SRscmos();
    thread_sr[num_dataset]->setFitParameters( min_photon, PSFSigma, zm_color, SRImage_pixelsize);
    thread_sr[num_dataset]->setThresholds( ll_threshold, subregion_cropx, subregion_cropy, N_photon_min, N_photon_max, loc_uncer_min, loc_uncer_max);

    int offset = 0;
    int frames_per_worker_fixed = frames_per_worker;
    if ( num_dataset < remaining_frames )
      frames_per_worker_fixed++;
    else
      offset = remaining_frames;

    int start = (frames_per_worker_fixed*num_dataset + offset); 
    int end = (frames_per_worker_fixed*num_dataset + offset + frames_per_worker_fixed ); 

    start += new_start;
    end += new_start;

    int total_num_frames = end - start;

    Image<uint16_t> *frame_stack;

#ifdef OMP
#pragma omp critical // can't read in parallel yet... 
#endif
    { 
      printf("(%2d) reading frame_range: %6d -- %6d [%d]\n",num_dataset, start, end, total_num_frames);
      frame_stack = input->parseImageData( start, end );
    }

    //printf("Analyzing %d files from %s\n", range, input->data_path); 
    thread_sr[num_dataset]->run(nm, frame_stack, num_dataset);
    if(frame_stack) delete frame_stack;
  }

#ifdef OMP
#pragma omp barrier
#endif

  SRscmos *sr = new SRscmos(*(thread_sr[0]));
  sr->setFitParameters( min_photon, PSFSigma, zm_color, SRImage_pixelsize);
  sr->setThresholds( ll_threshold, subregion_cropx, subregion_cropy, N_photon_min, N_photon_max, loc_uncer_min, loc_uncer_max);

  // merge results
  for ( int i = 0; i < number_of_iterations; i++ ) {
    SRscmos *one = thread_sr[i];

    int frames_per_worker_fixed = frames_per_worker;
    if ( i < remaining_frames )
      frames_per_worker_fixed++;

    sr->mapData( *(thread_sr[i]), i, frames_per_worker_fixed, remaining_frames);
    delete one;
  }

  if( createOutput( output_path ) != 0 ) {
    exit(EXIT_FAILURE);
  }

  ///// SAVE MATLAB FILE /////
  std::string output_mat_path = output_path;
  output_mat_path += "/data.mat";

  sr->writeMatFile( inversegain_path, varmap_path, data_path, output_mat_path.c_str(), input->number_of_frames );
  printf("saved data in %s file\n", output_mat_path.c_str());

  if ( sr->total_size == 0 ) {
    printf("nothing to plot\n");
    return;
  }

  std::string output_image_path = output_path;
  output_image_path += "out.tiff";

  sr->saveData( output_image_path.c_str() );
  printf("outputted reconstructed image \n");

  printf("plotting results (multithreaded)\n");
  plot( output_path.c_str(), sr); 

  delete nm;
  delete sr;

  delete input;
}

//SRscmos* ComputeManager::compute_subset(char* inversegain_path, char* varmap_path, char* data_path, char* output_path, 
//                                        int min_photon, float PSFSigma, int zm_color, float SRImage_pixelsize, int ll_threshold, 
//                                        int subregion_cropx, int subregion_cropy, int N_photon_min, int N_photon_max, float loc_uncer_min, float loc_uncer_max, 
//                                        int new_start, int new_end) {
//
//  /// Step 1: sCMOS camera noise /////
//  NoiseMap *nm = new NoiseMap();
//  nm->parseDataMat( inversegain_path, varmap_path );
//
//  ///// Step 2: Load Images //////
//  if ( !input )
//    input = new DataIO(data_path); // will know if it is a file or folder, and the number of frames
//  //DataIO *input = new DataIO(data_path); // will know if it is a file or folder, and the number of frames
//  input->number_of_frames = new_end - new_start - 1;
//  printf("Data contains: %d frames\n", input->number_of_frames);
//
//  int number_of_iterations = std::ceil((float)input->number_of_frames/MAX_FRAMES_PER_ITERATION);
//  number_of_iterations = std::max(MIN_THREAD_COUNT, number_of_iterations);
//
//  if ( input->number_of_frames < number_of_iterations )
//    number_of_iterations = input->number_of_frames;
//
//#ifdef OMP
//  omp_set_num_threads(number_of_iterations);
//  printf("Launching %d threads.\n", number_of_iterations);
//#endif
//
//  int frames_per_worker = floor( input->number_of_frames/number_of_iterations );
//  int remaining_frames = input->number_of_frames % number_of_iterations;
//
//  printf("Analyzing %d frames per worker.\n", frames_per_worker);
//  printf("Remaining %d frames, will be added to the first %d workers.\n", remaining_frames, remaining_frames);
//
//  SRscmos *sr = new SRscmos();
//  sr->setFitParameters( min_photon, PSFSigma, zm_color, SRImage_pixelsize);
//  sr->setThresholds( ll_threshold, subregion_cropx, subregion_cropy, N_photon_min, N_photon_max, loc_uncer_min, loc_uncer_max);
//
//  SRscmos *thread_sr[number_of_iterations];
//
//#ifdef OMP
//#pragma omp parallel for default(shared)
//#endif
//  for ( int num_dataset = 0; num_dataset < number_of_iterations; num_dataset++ ) {
//    //printf("\n------worker: %d-------\n", num_dataset);
//
//    thread_sr[num_dataset] = new SRscmos();
//    thread_sr[num_dataset]->setFitParameters( min_photon, PSFSigma, zm_color, SRImage_pixelsize);
//    thread_sr[num_dataset]->setThresholds( ll_threshold, subregion_cropx, subregion_cropy, N_photon_min, N_photon_max, loc_uncer_min, loc_uncer_max);
//
//    int offset = 0;
//    int frames_per_worker_fixed = frames_per_worker;
//    if ( num_dataset < remaining_frames )
//      frames_per_worker_fixed++;
//    else
//      offset = remaining_frames;
//
//    int start = (frames_per_worker_fixed*num_dataset + offset); 
//    int end = (frames_per_worker_fixed*num_dataset + offset + frames_per_worker_fixed ); 
//
//    start += new_start;
//    end += new_start;
//
//    int total_num_frames = end - start;
//
//    Frame256x256 *data_set;
//
//#ifdef OMP
//#pragma omp critical // can't read in parallel yet... 
//#endif
//    { 
//      printf("(%2d) reading frame_range: %6d -- %6d [%d]\n",num_dataset, start, end, total_num_frames);
//      data_set = (Frame256x256 *) malloc ( total_num_frames * sizeof(Frame256x256) );
//
//      ///// Step 0: Raw Data Acquisition /////
//      input->parseImageData( start, end, data_set );
//    }
//
//    //printf("Analyzing %d files from %s\n", range, input->data_path); 
//    thread_sr[num_dataset]->run(nm, data_set, num_dataset, total_num_frames);
//  }
//
//#ifdef OMP
//#pragma omp barrier
//#endif
//
//  // merge results
//  for ( int i = 0; i < number_of_iterations; i++ ) {
//    SRscmos *one = thread_sr[i];
//
//    int frames_per_worker_fixed = frames_per_worker;
//    if ( i < remaining_frames )
//      frames_per_worker_fixed++;
//
//    sr->mergeData( one->imtot, one->xtot, one->ytot, one->ztot, one->bgtot, one->lltot, one->photot, one->crlbxtot, one->crlbytot, one->maxprojim, one->sumim, one->total_size, i, frames_per_worker_fixed, remaining_frames );
//  }
//
//  return sr;
//
////  if( createOutput( output_path ) != 0 ) {
////    exit(EXIT_FAILURE);
////  }
////
////  ///// SAVE MATLAB FILE /////
////  std::string output_mat_path = output_path;
////  output_mat_path += "data.mat";
////
////  sr->writeMatFile( inversegain_path, varmap_path, data_path, output_path, input->number_of_frames );
////  printf("saved data in .mat file\n");
////
////  if ( sr->total_size == 0 ) {
////    printf("nothing to plot\n");
////    return;
////  }
////
////  std::string output_image_path = output_path;
////  output_image_path += "out.tiff";
////
////  sr->saveData( output_image_path.c_str() );
////  printf("outputted reconstructed image \n");
////
////  printf("plotting results (multithreaded)\n");
////  plot( output_path, sr); 
//}
//
//

void ComputeManager::compute(const std::string inversegain_path, const std::string varmap_path, const std::string data_path, const std::string output_path, 
                             int min_photon, float PSFSigma, int zm_color, float SRImage_pixelsize, int ll_threshold, 
                             int subregion_cropx, int subregion_cropy, int N_photon_min, int N_photon_max, float loc_uncer_min, float loc_uncer_max, 
                             int subset_size) {

  if ( !input ) 
    input = new DataIO(data_path); // will know if it is a file or folder, and the number of frames

  int num_of_iterations = input->number_of_frames / subset_size;

  for ( int i = 0; i <= num_of_iterations; i++ ) {
    int start = i*subset_size;
    int end = i*subset_size + subset_size; 

    if ( i == num_of_iterations ) {
      if ( input->number_of_frames % subset_size ) {
        start = i*subset_size;
        end = input->number_of_frames;
      }
    }

    //printf("%d -- %d\n", start, end);
 
    if( createOutput( output_path ) != 0 ) {
      exit(EXIT_FAILURE);
    }

    std::string output_path_subset = output_path;
    output_path_subset += "/" + std::to_string(start) + "_" + std::to_string(end) + "/";

    //printf("%s\n", output_path_subset.c_str());

    this->compute2( inversegain_path, varmap_path, data_path, (char*)output_path_subset.c_str(), min_photon, PSFSigma, zm_color, SRImage_pixelsize, ll_threshold, subregion_cropx, subregion_cropy, N_photon_min, N_photon_max, loc_uncer_min, loc_uncer_max, start, end);

  }
}

//// this is for 2-color data
//void ComputeManager::compute(char* inversegain_path, char* varmap_path, char* data_path, char* output_path, 
//                             int min_photon, float PSFSigma, int zm_color, float SRImage_pixelsize, int ll_threshold, 
//                             int subregion_cropx, int subregion_cropy, int N_photon_min, int N_photon_max, float loc_uncer_min, float loc_uncer_max, 
//                             int subset_size) {
//
//  if ( !input )
//    input = new DataIO(data_path); // will know if it is a file or folder, and the number of frames
//
//  int num_of_iterations = floor(input->number_of_frames / subset_size);
//
//  bool color = false; // false(0) == color1, true(1) == color2
//
//  int skip_frame_offset = 0;
//
//  int global_number_of_frames = input->number_of_frames;
//
//  SRscmos *sr1 = new SRscmos(); // master SR object.. will be merged here
//  sr1->setFitParameters( min_photon, PSFSigma, zm_color, SRImage_pixelsize);
//  sr1->setThresholds( ll_threshold, subregion_cropx, subregion_cropy, N_photon_min, N_photon_max, loc_uncer_min, loc_uncer_max);
//
//  SRscmos *sr2 = new SRscmos(); // master SR object.. will be merged here
//  sr2->setFitParameters( min_photon, PSFSigma, zm_color, SRImage_pixelsize);
//  sr2->setThresholds( ll_threshold, subregion_cropx, subregion_cropy, N_photon_min, N_photon_max, loc_uncer_min, loc_uncer_max);
//
//  SRscmos *thread_sr1[num_of_iterations/2]; // color1
//  int thread_sr1_count = 0;
//  SRscmos *thread_sr2[num_of_iterations/2]; // color2
//  int thread_sr2_count = 0;
//
//  SRscmos *sr_sub;
//
//  // compute super res for each color
//  for ( int i = 0; i <= num_of_iterations; i++ ) {
//    int start = i*subset_size + skip_frame_offset;
//    skip_frame_offset++;
//    int end = i*subset_size + subset_size + skip_frame_offset; 
//
//    if ( i == num_of_iterations ) {
//      //if ( input->number_of_frames % subset_size ) {
//        start = i*subset_size + skip_frame_offset;
//        end = global_number_of_frames;
//      //}
//   }
//
//    if ( start >= global_number_of_frames ) { // do not read non-existent frame(s)
//      break;
//    }
//
//    printf("(%d) %d -- %d\n", color, start, end);
//
//    if( createOutput( output_path ) != 0 ) {
//      exit(EXIT_FAILURE);
//    }
//
//    std::string output_path_subset = output_path;
//    output_path_subset += "/" + std::to_string(start) + "_" + std::to_string(end) + "/";
//
//    //printf("%s\n", output_path_subset.c_str());
//
//    sr_sub = this->compute_subset( inversegain_path, varmap_path, data_path, (char*)output_path_subset.c_str(), 
//                          min_photon, PSFSigma, zm_color, SRImage_pixelsize, ll_threshold,
//                          subregion_cropx, subregion_cropy, N_photon_min, N_photon_max, loc_uncer_min, loc_uncer_max, 
//                          start, end );
//
//    if ( color ) {
//      thread_sr1[thread_sr1_count++] = sr_sub;
//    }
//    else {
//      thread_sr2[thread_sr2_count++] = sr_sub; 
//    }
//
//    color = !color; // toggle color
//  }
//
//  // merge datasets for each color
//  for ( int i = 0; i < thread_sr1_count; i++ ) {
//    SRscmos *one = thread_sr1[i];
//
//    sr1->mergeData( one->imtot, one->xtot, one->ytot, one->ztot, one->bgtot, one->lltot, one->photot, one->crlbxtot, one->crlbytot, one->maxprojim, one->sumim, one->total_size, i, subset_size, 0 );
//  }
//
//  for ( int i = 0; i < thread_sr2_count; i++ ) {
//    SRscmos *one = thread_sr2[i];
//
//    sr2->mergeData( one->imtot, one->xtot, one->ytot, one->ztot, one->bgtot, one->lltot, one->photot, one->crlbxtot, one->crlbytot, one->maxprojim, one->sumim, one->total_size, i, subset_size, 0 );
//  }
//
//  if( createOutput( output_path ) != 0 ) {
//    exit(EXIT_FAILURE);
//  }
//
//  ///// SAVE MATLAB FILE /////
//  std::string output_mat_path = output_path;
//  output_mat_path += "/data_color1.mat";
//
//  sr1->writeMatFile( inversegain_path, varmap_path, data_path, output_mat_path.c_str(), input->number_of_frames );
//  printf("saved data in %s file\n", output_mat_path.c_str());
//
//  std::string output_image_path = output_path;
//  output_image_path += "sr_color1.tiff";
//
//  sr1->saveData( output_image_path.c_str() );
//  printf("outputted reconstructed image \n");
//
//
//  ///// SAVE MATLAB FILE /////
//  output_mat_path = output_path;
//  output_mat_path += "/data_color2.mat";
//
//  sr2->writeMatFile( inversegain_path, varmap_path, data_path, output_mat_path.c_str(), input->number_of_frames );
//  printf("saved data in %s file\n", output_mat_path.c_str());
//
//  output_image_path = output_path;
//  output_image_path += "sr_color2.tiff";
//
//  sr2->saveData( output_image_path.c_str() );
//  printf("outputted reconstructed image \n");
//}
