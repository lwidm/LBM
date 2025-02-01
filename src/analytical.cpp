#define _USE_MATH_DEFINES
#include "analytical.h"
#include <math.h>

static std::string filename = "analytical.cpp";

/**
 * \brief Provides the analytical solution to a 2D Poiseuille flow.
 *
 * This function calculates the state variables rho, ux and uy (not P) according
 * to the exact analytical solution to Poiseuille flow for a fixed value of the
 * external force `Fx`. Fx is obtained with the input parameters u_0, rho_0 and
 * nu. This function is used for setting initial conditions or estimating a
 * simulation's error.
 *
 * \param[in,out] state The state object that holds the flow variables rho, ux,
 * uy, and P (P not being modified here).
 * \param[in] gridsize The size of the computational grid in Nx Ny (Nz = 1, 2D
 * simulation).
 * \param[in] grid The grid object containing coordinate information.
 * \param[in] nu The kinematic viscosity of the fluid.
 * \param[in] rho_0 The initial density of the fluid.
 * \param[in] u_0 The characteristic velocity of the Poiseuille flow.
 * \param[in] p_0 The initial pressure inside the domain (a scalar).
 *
 * \see initial_condition
 */
void analytical_Poiseuille(State &state, const Gridsize &gridsize,
                           const Grid &grid, const double nu,
                           const double rho_0, const double u_0,
                           const double p_0) {
  const double Fx_val =
      12.0 * rho_0 * u_0 * nu / ((double)(gridsize[1]) * (double)(gridsize[1]));
  Eigen::ArrayXXd Fx =
      Eigen::ArrayXXd::Constant(gridsize[1], gridsize[0], Fx_val);
  state.rho = Eigen::ArrayXXd::Constant(gridsize[1], gridsize[0], rho_0);
  state.uy = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
  state.ux = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
  state.ux = Fx / (rho_0 * nu) * grid.Y * ((double)(gridsize[1]) - grid.Y);
  // BUG : Not sure if P is physical here
  state.P = Eigen::ArrayXXd::Constant(gridsize[1], gridsize[0], p_0);
}

/**
 * \brief Provides a 2D Taylor-Green flow initial_condition.
 *
 * This function calculates the state variables rho, ux and uy (not P) according
 * to the exact analytical solution of 2D Taylor-Green This function is used for
 * estimating a simulation's error.
 *
 * \param[in,out] state The state object that holds the flow variables rho, ux,
 * uy, and P (P not being modified here).
 * \param[in] gridsize The size of the computational grid in Nx Ny (Nz = 1, 2D
 * simulation).
 * \param[in] grid The grid object containing coordinate information.
 * \param[in] rho_0 The initial density of the fluid.
 * \param[in] u_0 The characteristic velocity of the Taylor-Green flow.
 * \param[in] p_0 The initial pressure inside the domain (a scalar).
 * solution
 *
 * \see initial_condition
 */
void initCond_TaylorGreen(State &state, const Gridsize &gridsize,
                          const Grid &grid, const double rho_0,
                          const double u_0, const double p_0) {
  const double L = std::min(gridsize[0], gridsize[1]);
  state.rho = Eigen::ArrayXXd::Constant(gridsize[1], gridsize[0], rho_0);
  state.ux =
      -u_0 * (2 * M_PI * grid.X / L).cos() * (2 * M_PI * grid.Y / L).sin();
  state.uy =
      u_0 * (2 * M_PI * grid.X / L).sin() * (2 * M_PI * grid.Y / L).cos();
  state.P =
      p_0 - (rho_0 * u_0 * u_0) / 4 *
                ((4 * M_PI * grid.X / L).cos() + (4 * M_PI * grid.Y / L).cos());
}

/**
 * \brief Provides the analytical solution to a 2D Taylor-Green flow.
 *
 * This function calculates the state variables rho, ux and uy (not P) according
 * to the exact analytical solution of 2D Taylor-Green This function is used for
 * estimating a simulation's error.
 *
 * \param[in,out] state The state object that holds the flow variables rho, ux,
 * uy, and P (P not being modified here).
 * \param[in] gridsize The size of the computational grid in Nx Ny (Nz = 1, 2D
 * simulation).
 * \param[in] grid The grid object containing coordinate information.
 * \param[in] nu The kinematic viscosity of the fluid.
 * \param[in] rho_0 The initial density of the fluid.
 * \param[in] u_0 The characteristic velocity of the Taylor-Green flow.
 * \param[in] p_0 The initial pressure inside the domain (a scalar).
 * \param[in] t_prime Dimensionless time at which to obtain the analytical
 * solution
 *
 * \see initial_condition
 */
void analytical_TaylorGreen(State &state, const Gridsize &gridsize,
                            const Grid &grid, const double nu,
                            const double rho_0, const double u_0,
                            const double p_0, const double t_prime) {
  const double L = std::min(gridsize[0], gridsize[1]);
  const double Re = u_0 * L / nu;
  initCond_TaylorGreen(state, gridsize, grid, rho_0, u_0, p_0);
  // BUG : Aliasing might be an issue here
  state.ux = state.ux * std::exp(-8 * M_PI * M_PI * t_prime / Re);
  state.uy = state.uy * std::exp(-8 * M_PI * M_PI * t_prime / Re);
  state.P = (state.P - p_0) * std::exp(-16 * M_PI * M_PI * t_prime / Re);
}
