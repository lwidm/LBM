#ifndef MAIN_H
#define MAIN_H

#include "Eigen/Dense"

typedef struct {
  Eigen::ArrayXXd rho;
  Eigen::ArrayXXd ux;
  Eigen::ArrayXXd uy;
  Eigen::ArrayXXd uz;
  Eigen::ArrayXXd P;
} State;

typedef std::array<std::size_t, 2> Gridsize;

typedef struct {
  Eigen::ArrayXd x;
  Eigen::ArrayXd y;
} GridVectors;

typedef struct {
  Eigen::ArrayXXd X;
  Eigen::ArrayXXd Y;
} Grid;

typedef enum { LBM, LBM_EXACT_DIFFERENCE, KBC } SolverType;

#endif // MAIN_H
