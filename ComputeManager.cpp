/**
*  @file    ComputeManager.cpp
*  @author  John Ravi (jjravi)
*  @date    5/13/2019
*  @version 2.3
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

  // default noise map
  inversegain_path = "./noise_map/la-gainmap.mat";
  varmap_path = "./noise_map/la-varmap.mat";

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
}

ComputeManager::~ComputeManager() {
  delete input;
}

void ComputeManager::compute(const std::string &data_path, const std::string &output_path, int new_start, int new_end) {

  /// Step 1: sCMOS camera noise ///
  NoiseMap *nm = new NoiseMap();
  nm->parseDataMat( inversegain_path, varmap_path );

  /// Step 2: Load Images ///
  if ( !input ) {
    input = new DataIO(data_path);
  }

  // create output folder if it does not exist
  if( createOutput( output_path ) != 0 ) {
    exit(EXIT_FAILURE);
  }

  for( int channel_index = 0; channel_index < input->getNumChannels(); channel_index++ ) {

    std::string channel_name = input->getChannelName(channel_index).c_str();
    printf("Analyzing channel(%d): %s\n", channel_index, channel_name.c_str());
    std::string channel_output_path = output_path + "/" + channel_name;

    // range of frames to analyze
    //new_end = (new_end == -1) ? input->number_of_frames : new_end;
    //new_start = (new_start == -1) ? 0 : new_start;
    //input->number_of_frames = new_end - new_start;

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
  //    start += new_start;
  //    end += new_start;

      int sub_num_frames = end - start;

      Image<uint16_t> *frame_stack;

#ifdef OMP
#pragma omp critical // faster to read sequentially
#endif
      {
        printf("(%2d) reading frame_range: %6d -- %6d [%d]\n",num_dataset, start, end, sub_num_frames);
        frame_stack = input->parseImageData( start, end, channel_index );
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

    if( createOutput( channel_output_path ) != 0 ) {
      exit(EXIT_FAILURE);
    }

    ///// SAVE MATLAB FILE /////
    std::string output_mat_path = channel_output_path;
    output_mat_path += "/" + input->file_name_noext + "_data.mat";
    sr->writeMatFile( inversegain_path, varmap_path, data_path, output_mat_path, input->number_of_frames );
    printf("saved data in %s file\n", output_mat_path.c_str());
    ////////////////////////////

    if ( sr->total_size == 0 ) {
      printf("nothing to plot\n");
      return;
    }

    ///// SR Reconstruction ////
    std::string output_image_path = channel_output_path;
    output_image_path += "/" + input->file_name_noext + "_sr.tiff";
    sr->saveData( output_image_path.c_str() );
    printf("outputted reconstructed image \n");
    ////////////////////////////

    printf("plotting results (multithreaded)\n");
    std::string plot_output_prefix = channel_output_path + "/" + input->file_name_noext;
    plot( plot_output_prefix.c_str(), sr);

    delete sr;

  }

  delete nm;
}

//void ComputeManager::compute(const std::string inversegain_path, const std::string varmap_path, const std::string data_path, const std::string output_path,
//                             int min_photon, float PSFSigma, int zm_color, float SRImage_pixelsize, int ll_threshold,
//                             int subregion_cropx, int subregion_cropy, int N_photon_min, int N_photon_max, float loc_uncer_min, float loc_uncer_max,
//                             int subset_size) {
//
//  if ( !input )
//    input = new DataIO(data_path); // will know if it is a file or folder, and the number of frames
//
//  //printf("number_of_frames: %d\n", input->number_of_frames);
//  //printf("subset_size: %d\n", subset_size);
//
//  int num_of_iterations = ceil(input->number_of_frames / subset_size);
//
//  for ( int i = 0; i < num_of_iterations; i++ ) {
//    int start = i*subset_size;
//    int end = i*subset_size + subset_size;
//
//    if ( i == num_of_iterations ) {
//      if ( input->number_of_frames % subset_size ) {
//        start = i*subset_size;
//        end = input->number_of_frames;
//      }
//    }
//
//    //printf("%d -- %d\n", start, end);
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
//    this->compute2( inversegain_path, varmap_path, data_path, (char*)output_path_subset.c_str(), min_photon, PSFSigma, zm_color, SRImage_pixelsize, ll_threshold, subregion_cropx, subregion_cropy, N_photon_min, N_photon_max, loc_uncer_min, loc_uncer_max, start, end);
//
//  }
//}

