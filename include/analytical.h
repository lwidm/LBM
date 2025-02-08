#ifndef ANALYTICAL_H
#define ANALYTICAL_H

#define _USE_MATH_DEFINES

#include "Eigen/Dense"
#include "lbm_core.h"

void analytical_Poiseuille_2D(State &state, const Gridsize &gridsize,
                              const Grid &grid, const double nu,
                              const double rho_0, const double u_0,
                              const double p_0);

void initCond_TaylorGreen_2D(State &state, const Gridsize &gridsize,
                             const Grid &grid, const double rho_0,
                             const double u_0, const double p_0);

void analytical_TaylorGreen_2D(State &state, const Gridsize &gridsize,
                               const Grid &grid, const double nu,
                               const double rho_0, const double u_0,
                               const double p_0, const double t_prime);

#endif // ANALYTICAL_H
