#ifndef HELPERS_H
#define HELPERS_H

#define __STDC_WANT_LIB_EXT1__ 1
#define __STDC_LIB_EXT1__
#include "Eigen/Dense"
#include "main.h"
#include <string>

Grid meshgrid(const Gridsize &gridsize, const GridVectors &gridvectors);

Eigen::ArrayXXd curlZ(const Eigen::ArrayXXd &ux, const Eigen::ArrayXXd &uy,
                      const Gridsize &gridsize, const double dr);

void print_err(std::string filename, std::string message);

#endif // HELPERS_H
