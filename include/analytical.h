#ifndef ANALYTICAL_H
#define ANALYTICAL_H

#include "Eigen/Dense"
#include "main.h"

void analytical_Poiseuille(State &state, const Gridsize &gridsize,
                           const Grid &grid, const double &nu,
                           const double &rho_0, const double &u_0);

#endif // ANALYTICAL_H
