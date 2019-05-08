/**
*  @file    main.cpp
*  @author  John Ravi (jjravi)
*  @date    5/7/2019
*  @version 1.3
*
*  @brief Provide CLI for run the FPALM computation.
*
*  @section DESCRIPTION
*
*  Initializes the computer manager and runs the simulation based on the given parameters.
*
*/

#include "ComputeManager.h"

int main(int argc, char* argv[]) {

  ComputeManager *cm = new ComputeManager();

  switch( argc ) {
    case 3: // compute using default parameters
      cm->compute(argv[1], argv[2]);
      break;

    case 16: // compute all frames
      cm->compute(argv[1], argv[2], argv[3], argv[4],
                   atoi(argv[5]), atof(argv[6]), atoi(argv[7]), atof(argv[8]),
                   atoi(argv[9]), atoi(argv[10]), atoi(argv[11]), atoi(argv[12]), atoi(argv[13]), atof(argv[14]), atof(argv[15]));
      break;

    case 17: // compute every <frame>
      cm->compute(argv[1], argv[2], argv[3], argv[4],
                   atoi(argv[5]), atof(argv[6]), atoi(argv[7]), atof(argv[8]),
                   atoi(argv[9]), atoi(argv[10]), atoi(argv[11]), atoi(argv[12]), atoi(argv[13]), atof(argv[14]), atof(argv[15]),
                   atoi(argv[16]));
      break;

    case 18: // compute <start_frame> <end_frame>
      cm->compute2(argv[1], argv[2], argv[3], argv[4],
                   atoi(argv[5]), atof(argv[6]), atoi(argv[7]), atof(argv[8]),
                   atoi(argv[9]), atoi(argv[10]), atoi(argv[11]), atoi(argv[12]), atoi(argv[13]), atof(argv[14]), atof(argv[15]),
                   atoi(argv[16]), atoi(argv[17]));
      break;

    default:
      printf("\n---Mode (default parameters)---\n");
      printf("usage: ./la-fpalm <data> <output-dir>\n");

      printf("\n---Mode ([start,end] frame)---\n");
      printf("usage: <program> <inversegain> <varmap> <data> <output-dir> <start_frame> <end_frame>\n");

      printf("\n---Mode (specify parameters)---\n");
      printf("usage: <program> <inversegain> <varmap> <data> <output-dir> "
              "<min_photon> <PSFSigma> <zm_color> <SRImage_pixelsize> "
              "<ll_threshold> <subregion_cropx> <subregion_cropy> <N_photon_min> <N_photon_max> <loc_uncer_min> <loc_uncer_max>\n");

      printf("\n---Mode (Batch every <subset_size> frames)---\n");
      printf("usage: <program> <inversegain> <varmap> <data> <output-dir> "
              "<min_photon> <PSFSigma> <zm_color> <SRImage_pixelsize> "
              "<ll_threshold> <subregion_cropx> <subregion_cropy> <N_photon_min> <N_photon_max> <loc_uncer_min> <loc_uncer_max>"
              "<subset_size> \n");

      return EXIT_FAILURE;
      break;
  }

  delete cm;

  printf("done\n");
  return EXIT_SUCCESS;
}
