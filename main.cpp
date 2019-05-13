/**
*  @file    main.cpp
*  @author  John Ravi (jjravi)
*  @date    5/13/2019
*  @version 2.0
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

  if( argc < 2 ) {
    printf("usage: <program> <mode> ...\n");
    printf("\n--- Mode 0, analyze all frames ---\n");
    //printf("\n--- Mode 1, analyze [start, end] frames ---\n");

    return EXIT_FAILURE;
  }

  switch( atoi(argv[1]) ) {
    case 0: // compute using default parameters
      printf("\n--- Mode 0, analyze all frames ---\n");

      if( argc == 4 ) {
        cm->compute(argv[2], argv[3]);
      }
      else if( argc == 17 ) {
        cm->inversegain_path = argv[2];
        cm->varmap_path = argv[3];
 
        cm->min_photon = atoi(argv[6]);
        cm->PSFSigma = atof(argv[7]);
        cm->zm_color = atoi(argv[8]);
        cm->SRImage_pixelsize = atof(argv[9]);
      
        cm->ll_threshold = atoi(argv[10]); 
        cm->subregion_cropx = atoi(argv[11]); 
        cm->subregion_cropy = atoi(argv[12]);
        cm->N_photon_min = atoi(argv[13]);
        cm->N_photon_max = atoi(argv[14]);
        cm->loc_uncer_min = atof(argv[15]);
        cm->loc_uncer_max = atof(argv[16]);

        cm->compute(argv[4], argv[5]);
      }
      else {
        // (default parameters)
        printf("usage: <program> 0 <data> <output-dir>\n");

        // (specify parameters)
        printf("usage: <program> 0 <inversegain> <varmap> <data> <output-dir> "
                "<min_photon> <PSFSigma> <zm_color> <SRImage_pixelsize> "
                "<ll_threshold> <subregion_cropx> <subregion_cropy> <N_photon_min> <N_photon_max> <loc_uncer_min> <loc_uncer_max>\n");
 
        return EXIT_FAILURE;
      }

      break;

//    case 1: // compute <start_frame> <end_frame>
//
//      if( argc !=  
//        printf("\n--- Mode 1, analyze [start, end] frames ---\n");
//        printf("usage: <program> <inversegain> <varmap> <data> <output-dir> <start_frame> <end_frame>\n");
//        return EXIT_FAILURE;
//      }
//
//      cm->compute2(argv[1], argv[2], argv[3], argv[4],
//                   atoi(argv[5]), atof(argv[6]), atoi(argv[7]), atof(argv[8]),
//                   atoi(argv[9]), atoi(argv[10]), atoi(argv[11]), atoi(argv[12]), atoi(argv[13]), atof(argv[14]), atof(argv[15]),
//                   atoi(argv[16]), atoi(argv[17]));
//      break;
//
//    case 2: // compute every <frame>
//      if( argc != 
//        printf("\n--- Mode 2, Batch every <subset_size> frames ---\n");
//        printf("usage: <program> <inversegain> <varmap> <data> <output-dir> "
//                "<min_photon> <PSFSigma> <zm_color> <SRImage_pixelsize> "
//                "<ll_threshold> <subregion_cropx> <subregion_cropy> <N_photon_min> <N_photon_max> <loc_uncer_min> <loc_uncer_max>"
//                "<subset_size> \n");
//        return EXIT_FAILURE;
//      }
//
//      cm->compute(argv[1], argv[2], argv[3], argv[4],
//                   atoi(argv[5]), atof(argv[6]), atoi(argv[7]), atof(argv[8]),
//                   atoi(argv[9]), atoi(argv[10]), atoi(argv[11]), atoi(argv[12]), atoi(argv[13]), atof(argv[14]), atof(argv[15]),
//                   atoi(argv[16]));
//      break;

    default:
      printf("\n---Unsupported Mode---\n");
      return EXIT_FAILURE;

      break;
  }

  delete cm;

  printf("done\n");
  return EXIT_SUCCESS;
}
