#include "main.h"
#include "Eigen/src/Core/Array.h"
#include "analytical.h"
#include "export_data.h"
#include "helpers.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <sstream>

static std::string this_filename = "main.cpp";
const unsigned int D = 2;
const unsigned int Q = 9;

/**
 * \brief Initializes the state variables for the lattice Boltzmann simulation.
 *
 * This function sets the initial conditions for the lattice Boltzmann
 * simulation. The user needs to adjust the internal function call to the
 * desired initial condition function, such as `analytical_Poiseuille`, or any
 * other suitable function based on their specific requirements.
 *
 * \param[in,out] state The state object that holds the flow variables rho, ux,
 * uy, and P.
 * \param[in] gridsize The size of the computational grid in Nx Ny (Nz = 1, 2D
 * simulation).
 * \param[in] grid The grid object containing coordinate information (grid.X,
 * grid.Y <- 2D eigen arrays).
 * \param[in] nu The kinematic viscosity of the fluid.
 * \param[in] rho_0 The initial density of the fluid.
 * \param[in] u_0 The characteristic velocity of the flow.
 * \param[in] p_0 The initial pressure inside the domain (a scalar).
 */
void initial_condition(State &state, const Gridsize &gridsize, const Grid &grid,
                       const double &nu, const double &rho_0, const double &u_0,
                       const double &p_0) {
  // analytical_Poiseuille(state, gridsize, grid, nu, rho_0, u_0, p_0);
  initCond_TaylorGreen(state, gridsize, grid, nu, rho_0, u_0, p_0);
}

/**
 * \brief The streaming step of the Lattice Baltzmann method
 *
 * This function applies the streaming to the populations f. To determine in
 * which direction every population f[i] should be streamed the vecotors csx and
 * csy are used. f cxs and cys are array of length Q where Q is the number of
 * populations in the simulation (e.g. D2Q9 -> Q = 9, D1Q3 -> Q = 3). The
 * cirshift function is used for this computation.
 *
 * \param[in, out] f The populations to be shifted
 * \param[in] gridsize The size of the grid in the format [Nx, Ny, Nz] (Nz = 1)
 * \param[in] cxs The x component in which each population f[i] should be
 * streamed
 * \param[in] cys The y component in which each population f[i] should be
 * streamed
 *
 * \see circshift: Cirshift operation on Eigen Arrays
 */
void streaming_step(std::array<Eigen::ArrayXXd, Q> &f, const Gridsize &gridsize,
                    const std::array<double, Q> &cxs,
                    const std::array<double, Q> &cys) {
  for (unsigned int i = 1; i < Q; ++i) {
    roll2D(f[i], gridsize, cys[i], cxs[i]);
  }
}

/**
 * \brief Compute macroscopic variables (e.g. rho, ux, uy, uz, P) from the
 * populations f
 *
 * This function computes the macroscopic variables rho, ux, uy and uz from the
 * populations f
 *
 * \param[in, out] state The state of the simulations which is a struct that
 * hold the arrays for rho, ux, uy, uz and P (i.e. the state)
 * \param[in] f The current population distribution of the simulation
 * \param[in] gridsize The size of the grid in the format [Nx, Ny, Nz] (Nz = 1)
 * \param[in] cxs The x component of the grid velocities associated with each
 * f[i]
 * \param[in] cys The y component of the grid velocities associated with each
 * f[i]
 * \param[in] cs The lattice speed of sound
 */
void compute_macroscopic_variables(State &state, const Gridsize &gridsize,
                                   const std::array<Eigen::ArrayXXd, Q> &f,
                                   const std::array<double, Q> &cxs,
                                   const std::array<double, Q> &cys,
                                   const double &cs) {
  state.rho = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
  state.ux = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
  state.uy = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
  state.P = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
  for (unsigned int i = 0; i < Q; ++i) {
    state.rho += f[i];
    state.ux += cxs[i] * f[i];
    state.uy += cys[i] * f[i];
  }
  state.ux = state.ux / state.rho;
  state.uy = state.uy / state.rho;
  // BUG : Not sure about pressure term
  state.P = cs * cs + (state.ux * state.ux + state.uy * state.uy);
}

/**
 * \brief Main simulation function for the lattice Boltzmann method.
 *
 * This function runs the lattice Boltzmann simulation for the specified
 * simulation time. It modifies the state variables rho, ux, (uy, uz) and P
 * according to the chosen solver type. The function initializes the state using
 * the provided initial condition function and then proceeds with the
 * simulation.
 *
 * \param[in,out] state The state object that holds the flow variables rho, ux,
 * (uy, uz) and P.
 * \param[in] gridsize The size of the computational grid in Nx, Ny, Nz (Nz = 1
 * for 2D simulations).
 * \param[in] grid The grid object containing coordinate information (1D/2D/3D
 * arrays of X, X,Y, X,Y,Z depending on the simulation dimension).
 * \param[in] nu The viscosity of the fluid
 * \param[in] sim_time The total simulation time for which the solver should
 * run.
 * \param[in] solver An enum representing the type of lattice Boltzmann solver
 * to be used (e.g., LBM, LBM_EXACT_DIFFERENCE, or KBC). solver to be used
 * (e.g., "LBM", "KBC", or "LBM_exact_difference").
 * \param[in] initial_condition A function that sets the initial conditions for
 * the simulation.
 * \return true if the entire simulation was run for the specified sim_time,
 * false if it terminated early.
 *
 * \details The function returns false if the simulation terminates early. Early
 * termination are implementation specific and depend on the users
 * implementation. For e.g. one may want the function to terminate early if the
 * simulation becomes unstable or convergences.
 */
int lattice_bolzmann_simulation(
    State &state, const Gridsize &gridsize, const Grid &grid, const double &nu,
    const unsigned int &sim_time, const SolverType &solver,
    std::function<void(State &, const Gridsize &, const Grid &)>
        initial_condition) {

  // output definitions
  const unsigned int output_freq = 20;
  const bool save_output = true;
  int err = 0;

  // D2Q9 simulation
  const unsigned int dr = 1; // BUG : Not sure if this should depend on grid
  const unsigned int dt = 1; // BUG : Not sure if this should be an input
  const double c = (double)dr / (double)dt;
  const double cs = c / std::sqrt(3.);
  const double beta = (double)dt / (2 * nu / (cs * cs) + (double)dt);

  const std::array<double, Q> cxs = {0, 0, c, c, c, 0, -c, -c, -c};
  const std::array<double, Q> cys = {0, c, c, 0, -c, -c, -c, 0, c};
  const std::array<double, Q> weigths = {4. / 9,  1. / 9,  1. / 36,
                                         1. / 9,  1. / 36, 1. / 9,
                                         1. / 36, 1. / 9,  1. / 36};
  auto feq_calc = [&cxs, &cys, &weigths, &c,
                   &cs](std::array<Eigen::ArrayXXd, Q> &feq,
                        const Eigen::ArrayXXd &rho, const Eigen::ArrayXXd &ux,
                        const Eigen::ArrayXXd &uy) {
    for (unsigned int i = 0; i < Q; ++i) {
      feq[i] = rho * weigths[i] *
               (1.0 + 3. * (cxs[i] * ux + cys[i] * uy) +
                9. * (cxs[i] * ux + cys[i] * uy).square() / 2 -
                3 * (ux.square() + uy.square()) / 2);
    }
  };

  // Initial condition
  initial_condition(state, gridsize, grid);
  std::array<Eigen::ArrayXXd, Q> f;
  std::array<Eigen::ArrayXXd, Q> feq;
  std::array<Eigen::ArrayXXd, Q> fmirr;
  feq_calc(feq, state.rho, state.ux, state.uy);
  f = feq;

  if (save_output) {
    if (save_state("taylor green", "num_", state, gridsize, grid, (double)0,
                   CONFIRM) != 0) {
      LOG_ERR(this_filename, "Saving state failed");
      return 1;
    }
  }

  for (unsigned int t = 1; t <= sim_time; ++t) {
    // Collision step
    feq_calc(feq, state.rho, state.ux, state.uy);
    for (unsigned int i = 0; i < Q; ++i) {
      fmirr[i] = 2 * feq[i] - f[i];
      f[i] = (1 - beta) * f[i] + beta * fmirr[i];
    }
    // Streaming step
    streaming_step(f, gridsize, cxs, cys);
    // Boundary Condition
    // Compute macroscopic variables
    compute_macroscopic_variables(state, gridsize, f, cxs, cys, cs);

    if (save_output && t % output_freq == 0) {
      if (save_state("taylor green", "num_", state, gridsize, grid, (double)t,
                     CONFIRM) != 0) {
        LOG_ERR(this_filename, "Saving state failed");
        return 1;
      }
    }
    std::cout << "\r time: " << t << "/" << sim_time;
  }

  return 0;
}

/**
 * \brief Sets up and runs the lattice Boltzmann simulation.
 *
 * This function is the entry point of the program. It allows the user to set up
 * the simulation parameters according to their specific needs, run the
 * simulation and save/postprocess the results.
 *
 * \return 0 if the program completed without (critical) errors, a non-zero
 * value otherwise.
 *
 */
int main() {
  // ----------- flow setup -----------
  const double dr = 1;
  const double dt = 1;
  // const unsigned int sim_time = 1000;
  const double rho_0 = 1;
  const double p_0 = rho_0 / 3;
  const double Re = 100;
  // const double nu = 0.064;
  const double u_0 = 0.1;
  const double t_prime = 2.5;
  const double L = 128;
  const double nu = u_0 * L / Re;
  const unsigned int sim_time = (unsigned int)(t_prime * L / u_0);

  // // ----------- Run the 3 simulations -----------

  MetaData metadata;
  std::ostringstream oss;
  oss << "taylor green";
  std::string sim_name = oss.str();
  int err;
  err = init_save_dir(sim_name, metadata, CONFIRM);
  if (err != 0) {
    LOG_ERR(this_filename, "Createing new save failed");
  }

  // for (std::size_t i = 0; i < 3; ++i) {

  // std::size_t Nx = 22 * D[i];
  // std::size_t Ny = (std::size_t)(4.1 * (double)(D[i]));
  std::size_t Nx = L;
  std::size_t Ny = L;
  std::size_t Nz = 1;

  const Gridsize gridsize = std::array<std::size_t, 3>{Nx, Ny, Nz};

  GridVectors gridvectors;
  gridvectors.x = Eigen::ArrayXd::LinSpaced(Nx, dr / 2, Nx * dr);
  gridvectors.y = Eigen::ArrayXd::LinSpaced(Ny, dr / 2, Ny * dr);
  Grid grid = meshgrid(gridsize, gridvectors);

  auto initial_condition_preset = [&nu, &rho_0, &u_0,
                                   &p_0](State &state, const Gridsize &gridsize,
                                         const Grid &grid) {
    initial_condition(state, gridsize, grid, nu, rho_0, u_0, p_0);
  };

  State state;
  int stable = lattice_bolzmann_simulation(state, gridsize, grid, nu, sim_time,
                                           LBM, initial_condition_preset);

  Eigen::ArrayXXd curl = curlZ(state.ux, state.uy, gridsize, dr);

  State analytical_state;
  analytical_TaylorGreen(analytical_state, gridsize, grid, nu, rho_0, u_0, p_0,
                         t_prime);
  if (save_state("taylor green", "ana_", analytical_state, gridsize, grid,
                 (double)sim_time, CONFIRM) != 0) {
    LOG_ERR(this_filename, "Saving state failed");
    return 1;
  }

  // }

  return 0;
}
