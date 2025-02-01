#include "main.h"
#include "analytical.h"
#include "export_data.h"
#include "helpers.h"
#include "lbm_core.h"

#include <array>
#include <cstddef>
#ifdef GOOGLE_PROF
#include <gperftools/profiler.h>
#endif
#include <iomanip>
#include <sstream>
#include <string>

static std::string this_filename = "main.cpp";

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
                       [[maybe_unused]] const double nu, const double rho_0,
                       const double u_0, const double p_0) {
  // analytical_Poiseuille(state, gridsize, grid, nu, rho_0, u_0, p_0);
  initCond_TaylorGreen(state, gridsize, grid, rho_0, u_0, p_0);
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
#ifdef GOOGLE_PROF
  std::cout << "using goolge profiler" << std::endl;
  ProfilerStart("LBM.prof");
#endif
  // ----------- flow setup -----------
  const double dr = 1;
  // const double dt = 1;
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
  const SolverType solver = LBM;

  // // ----------- Run the 3 simulations -----------

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

  MetaData metadata;
  std::string sim_name = "2D_Taylor_Green";
  std::string description =
      "A standard LBM simulation of a 2D Taylor Green flow. Since the "
      "analytical solution of a Taylor Green flow is know this can be used to "
      "calculate the simulations' error.";
  std::string initial_condition_string = "2D Taylor Green";
  create_metadata(metadata, 1, sim_name, description, gridsize, nu, rho_0, u_0,
                  p_0, solver, initial_condition_string);
  if (init_save_dir(sim_name, metadata, grid, CONFIRM) != 0) {
    LOG_ERR(this_filename, "Creating new save failed");
  }

  auto initial_condition_preset = [&nu, &rho_0, &u_0,
                                   &p_0](State &state, const Gridsize &gridsize,
                                         const Grid &grid) {
    initial_condition(state, gridsize, grid, nu, rho_0, u_0, p_0);
  };

  State state;
  int stable =
      lattice_bolzmann_simulation(state, sim_name, gridsize, grid, nu, sim_time,
                                  solver, initial_condition_preset);

  Eigen::ArrayXXd curl = curlZ(state.ux, state.uy, gridsize, dr);

  State analytical_state;
  analytical_TaylorGreen(analytical_state, gridsize, grid, nu, rho_0, u_0, p_0,
                         t_prime);
  if (save_state(sim_name, "ana_", analytical_state, (double)sim_time,
                 CONFIRM) != 0) {
    LOG_ERR(this_filename, "Saving state failed");
    return 1;
  }

  // }

#ifdef GOOGLE_PROF
  std::cout << "using goolge profiler" << std::endl;
  ProfilerStop();
#endif
  return stable;
}
