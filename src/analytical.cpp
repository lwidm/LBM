#include "analytical.h"

static std::string filename = "analytical.cpp";

/**
 * \brief Provides the analytical solution to a 2D Poiseuille flow.
 *
 * This function calculates the state variables rho, ux and uy (not P) according
 * to the exact analytical solution to Poiseuille flow for a fixed value of the
 * external force `Fx`, calculated based on the input parameters. It is used for
 * setting initial conditions or estimating the simulation's error.
 *
 * \param[in,out] state The state object that holds the flow variables rho, ux,
 * uy, and P (P not being modified here).
 * \param[in] gridsize The size of the computational grid in Nx Ny (Nz = 1, 2D
 * simulation).
 * \param[in] grid The grid object containing coordinate information.
 * \param[in] nu The kinematic viscosity of the fluid.
 * \param[in] rho_0 The initial density of the fluid.
 * \param[in] u_0 The characteristic velocity of the Poiseuille flow.
 */
void analytical_Poiseuille(State &state, const Gridsize &gridsize,
                           const Grid &grid, const double &nu,
                           const double &rho_0, const double &u_0) {
  const double Fx_val =
      12.0 * rho_0 * u_0 * nu / ((double)(gridsize[1]) * (double)(gridsize[1]));
  Eigen::ArrayXXd Fx =
      Eigen::ArrayXXd::Constant(gridsize[1], gridsize[0], Fx_val);
  state.rho = Eigen::ArrayXXd::Constant(gridsize[1], gridsize[0], rho_0);
  state.uy = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
  state.ux = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
  state.ux = Fx / (rho_0 * nu) * grid.Y * ((double)(gridsize[1]) - grid.Y);
}
