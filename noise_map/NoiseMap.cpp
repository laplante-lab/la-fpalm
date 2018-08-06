#include "NoiseMap.h"

NoiseMap::NoiseMap() {
  gain_map = NULL;
  var_map = NULL;
  o_map = NULL;
}

NoiseMap::~NoiseMap() {
  if ( gain_map ) { 
    delete gain_map;
  }

  if ( var_map ) {
    delete var_map;
  }

  if ( o_map ) {
    delete o_map;
  }
}

void NoiseMap::parseDataMat( const std::string &gain_path, const std::string &var_path ) {

  mat_t *matgainfp = NULL;
  mat_t *matvarfp = NULL;
  matvar_t *matgainvar = NULL;
  matvar_t *matovar = NULL;
  matvar_t *matvarvar = NULL;

  matgainfp = Mat_Open(gain_path.c_str(),MAT_ACC_RDONLY);
  matvarfp = Mat_Open(var_path.c_str(),MAT_ACC_RDONLY);

  if ( matgainfp == NULL ) {
    fprintf(stderr,"Error opening MAT file %s\n", gain_path.c_str());
    exit( EXIT_FAILURE );
  }

  if ( matvarfp == NULL ) {
    fprintf(stderr,"Error opening MAT file %s\n", var_path.c_str());
    exit( EXIT_FAILURE );
  }

  matgainvar = Mat_VarReadInfo(matgainfp,"inversegain");
  matovar = Mat_VarReadInfo(matvarfp,"omap");
  matvarvar = Mat_VarReadInfo(matvarfp,"varmap");

  if ( matgainvar == NULL ) {
    printf("inversegain was not found in %s\n", gain_path.c_str());
  }

  if ( matovar == NULL || matvarvar == NULL ) {
    printf("omap or varmap was not found in %s\n", var_path.c_str());
  }

  //Mat_VarPrint(matgainvar, 1);
  //Mat_VarPrint(matovar, 1);
  //Mat_VarPrint(matvarvar, 1);

  int start[2]={0,0};
  int stride[2]={1,1};
  int edge[2];

  if ( !gain_map ) {
    gain_map = new Image<double>(matgainvar->dims[0], matgainvar->dims[1]);
  }

  if ( !var_map ) {
    var_map = new Image<double>(matovar->dims[0], matovar->dims[1]);
  }

  if ( !o_map ) {
    o_map = new Image<double>(matvarvar->dims[0], matvarvar->dims[1]);
  }

  edge[0] = matgainvar->dims[0];
  edge[1] = matgainvar->dims[1];

  Mat_VarReadData( matgainfp, matgainvar, gain_map->pixels, start, stride, edge);

  edge[0] = matovar->dims[0];
  edge[1] = matovar->dims[1];
 
  Mat_VarReadData( matvarfp, matovar, o_map->pixels, start, stride, edge);

  edge[0] = matvarvar->dims[0];
  edge[1] = matvarvar->dims[1];
 
  Mat_VarReadData( matvarfp, matvarvar, var_map->pixels, start, stride, edge);

//  // crop out subsection of camera noise
//  int mapbase = IMAGE_SIZE - IMAGE_SIZE/2;
//
//  for( int i = mapbase, indexx = 0; i < mapbase+IMAGE_SIZE; i++, indexx++ ) {
//    for( int j = mapbase, indexy = 0; j < mapbase+IMAGE_SIZE; j++, indexy++ ) {
//    //  printf("%0.4lf ", full_gain_map[i*512 + j]);
//
//      gain_map[indexx*256 + indexy] = full_gain_map[ j*512 + i ];
//      var_map[indexx*256 + indexy]  = full_var_map[ j*512 + i ];
//      o_map[indexx*256 + indexy]    = full_o_map[ j*512 + i ];
//
//    }
////    printf("%0.4lf ", full_gain_map[0*512 + i]);
////    printf("\n");
//  }
//
//  //printf("\n\n\n");
//
//  //for ( int i = 0; i < 256; i++ ) {
//  //  printf("i(%d):%0.4lf\n", i, gain_map[0*256 + i]);
//
//  //}

  Mat_VarFree(matgainvar);
  Mat_VarFree(matovar);
  Mat_VarFree(matvarvar);

  Mat_Close(matgainfp);
  Mat_Close(matvarfp);
}

// Deprecated (0.5s speedup but need to convert to binary file)
//void NoiseMap::parseData( char *gain_path, char *var_path, char *omap_path ) {
//  FILE * gain_cal_file = fopen(gain_path,"rb");
//  FILE * var_cal_file = fopen(var_path,"rb");
//  FILE * omap_file = fopen(omap_path,"rb");
//
//  if ( gain_cal_file==NULL || var_cal_file==NULL || omap_file==NULL ) {
//  printf("cannot open noise_maps\n");
//  exit(1);
//  }
//
//  // CHECK THE SIZES HERE!
//  double *full_gain_map = (double *)malloc( NOISE_MAP_SIZE * NOISE_MAP_SIZE * sizeof(double) );
//  double *full_var_map = (double *)malloc( NOISE_MAP_SIZE * NOISE_MAP_SIZE * sizeof(double) );
//  double *full_o_map = (double *)malloc( NOISE_MAP_SIZE * NOISE_MAP_SIZE * sizeof(double) );
//
//  fread( full_gain_map, sizeof( double ), NOISE_MAP_SIZE*NOISE_MAP_SIZE, gain_cal_file );
//  fread( full_var_map, sizeof( double ), NOISE_MAP_SIZE*NOISE_MAP_SIZE, var_cal_file );
//  fread( full_o_map, sizeof( double ), NOISE_MAP_SIZE*NOISE_MAP_SIZE, omap_file );
//
//  fclose (gain_cal_file);
//  fclose (var_cal_file);
//  fclose (omap_file);
//
//  // crop out subsection of camera noise
//  int mapbase = IMAGE_SIZE - IMAGE_SIZE/2;
//
//  for ( int i = mapbase, index = 0; i < mapbase+IMAGE_SIZE; i++, index++ ) {
//  std::memcpy((gain_map+index*IMAGE_SIZE), (full_gain_map+i*NOISE_MAP_SIZE+mapbase), IMAGE_SIZE*sizeof(double));
//  std::memcpy((var_map+index*IMAGE_SIZE), (full_var_map+i*NOISE_MAP_SIZE+mapbase), IMAGE_SIZE*sizeof(double));
//  std::memcpy((o_map+index*IMAGE_SIZE), (full_o_map+i*NOISE_MAP_SIZE+mapbase), IMAGE_SIZE*sizeof(double));
//  }
//
//  free( full_gain_map );
//  free( full_var_map );
//  free( full_o_map );
//}

