#ifndef HELPERS_H
#define HELPERS_H

#define __STDC_WANT_LIB_EXT1__ 1
#define __STDC_LIB_EXT1__
#include "Eigen/Dense"
#include "main.h"

// clang-format off
typedef enum {
  LOG_LEVEL_OFF = 0,
  LOG_LEVEL_ERR,
  LOG_LEVEL_WRN,
  LOG_LEVEL_INF,
  LOG_LEVEL_DBG
} Log_level;
// clang-format on

void log(std::string filename, std::string message, Log_level log_level);

#define LOG_LEVEL LOG_LEVEL_DBG
#define LOG_ERR(filename, message) log(filename, message, LOG_LEVEL_ERR)
#define LOG_WRN(filename, message) log(filename, message, LOG_LEVEL_WRN)
#define LOG_INF(filename, message) log(filename, message, LOG_LEVEL_INF)
#define LOG_DBG(filename, message) log(filename, message, LOG_LEVEL_DBG)

Grid meshgrid(const Gridsize &gridsize, const GridVectors &gridvectors);

Eigen::ArrayXXd curlZ(const Eigen::ArrayXXd &ux, const Eigen::ArrayXXd &uy,
                      const Gridsize &gridsize, const double dr);

#endif // HELPERS_H
