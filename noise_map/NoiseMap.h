#ifndef NOISE_MAP_H
#define NOISE_MAP_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "matio.h"

#include "../parameters.h"

#include "../Image.h"

class NoiseMap {

  public:
    Image<double> *gain_map;
    Image<double> *var_map;
    Image<double> *o_map;

    NoiseMap();
    ~NoiseMap(); 
    void parseDataMat( const std::string &gain_path, const std::string &var_path );
};

#endif
