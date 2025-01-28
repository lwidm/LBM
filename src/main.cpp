#include "main.h"
#include "analytical.h"
#include "export_data.h"
#include "helpers.h"

#include <array>
#include <cstddef>
#include <functional>
#include <sstream>

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
 */
void initial_condition(State &state, const Gridsize &gridsize, const Grid &grid,
                       const double &nu, const double &rho_0,
                       const double &u_0) {
  analytical_Poiseuille(state, gridsize, grid, nu, rho_0, nu);
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
int lattice_bolzmann_solver(
    State &state, const Gridsize &gridsize, const Grid &grid,
    const unsigned int &sim_time, const SolverType &solver,
    std::function<void(State &, const Gridsize &, const Grid &)>
        initial_condition) {
  initial_condition(state, gridsize, grid);
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
  double dr = 1;
  double dt = 1;
  double sim_time = 1000;
  double rho_0 = 1;
  double p_0 = rho_0 / 3;
  double Re = 20;
  double nu = 0.02;
  std::array<std::size_t, 3> D = {12, 16, 20};
  std::array<double, 3> u_0;
  for (std::size_t i = 0; i < 3; ++i) {
    u_0[i] = Re / (double)(D[i]) * nu;
  }

  // ----------- Run the 3 simulations -----------
  for (std::size_t i = 0; i < 3; ++i) {

    std::size_t Nx = 22 * D[i];
    std::size_t Ny = (std::size_t)(4.1 * (double)(D[i]));
    std::size_t Nz = 1;

    const Gridsize gridsize = std::array<std::size_t, 3>{Nx, Ny, Nz};

    GridVectors gridvectors;
    gridvectors.x = Eigen::ArrayXd::LinSpaced(Nx, dr / 2, Nx * dr);
    gridvectors.y = Eigen::ArrayXd::LinSpaced(Ny, dr / 2, Ny * dr);
    Grid grid = meshgrid(gridsize, gridvectors);

    auto initial_condition_preset = [&nu, &rho_0, &u_0,
                                     i](State &state, const Gridsize &gridsize,
                                        const Grid &grid) {
      initial_condition(state, gridsize, grid, nu, rho_0, u_0[i]);
    };

    State state;
    int stable = lattice_bolzmann_solver(state, gridsize, grid, sim_time, LBM,
                                         initial_condition_preset);

    Eigen::ArrayXXd curl = curlZ(state.ux, state.uy, gridsize, dr);

    MetaData metadata;
    std::ostringstream oss;
    oss << "test" << i;
    std::string sim_name = oss.str();
    int err;
    err = init_save_dir(sim_name, metadata, CONFIRM);
    if (err != 0) {
      LOG_ERR(this_filename, "Createing new save failed");
    }
    err = save_state(sim_name, state, gridsize, grid, sim_time, CONFIRM);
    if (err != 0) {
      LOG_ERR(this_filename, "Saving state failed");
    }
  }

  LOG_ERR(this_filename, "Test ERR");
  LOG_WRN(this_filename, "Test WRN");
  LOG_INF(this_filename, "Test INF");
  LOG_DBG(this_filename, "Test DBG");

  return 0;
}
