#include "main.h"
#include "analytical.h"
#include "helpers.h"
#include "py_action.h"

#include <array>
#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>

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
 * \details The function returns false if the simulation terminates early due to
 * reasons such as instability, convergence, or other case-specific user set
 * conditions.
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
 * \return 0 if the simulation runs successfully, a non-zero value otherwise.
 *
 */
int main() {
  // ----------- flow setup -----------
  double dr = 1;
  double dt = 1;
  double rho_0 = 1;
  double p_0 = rho_0 / 3;
  double Re = 20;
  double nu = 0.02;
  std::array<std::size_t, 3> D = {12, 16, 20};
  std::array<double, 3> u_0;
  for (std::size_t i = 0; i < 3; ++i) {
    u_0[i] = Re / (double)(D[i]) * nu;
  }

  // ----------- initializing shared memory for python actions -----------
  auto p_sharedmemory_struct =
      unique_ptr_shared_memory(new SharedMemory_struct{nullptr, nullptr});
  std::shared_ptr<double[]> shared_dArray;
  std::shared_ptr<int> semph_ptr;
  FILE *pipe;
  std::size_t max_size = (std::size_t)((22.0 * (double)(D[2])) *
                                       (std::size_t)(4.1 * (double)(D[2]))) +
                         1;
  std::size_t mmap_max_vector_size = max_size + max_size + max_size;
  int err = py_init(p_sharedmemory_struct, shared_dArray, semph_ptr, pipe,
                    mmap_max_vector_size);
  std::cout << "exit of py_init" << std::endl;
  if (err != 0) {
    print_err("main.cpp", "Failed to initialize py_action");
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
    // int stable = lattice_bolzmann_solver(state, gridsize, grid, LBM,
    //                                      initial_condition_preset);

    initial_condition(state, gridsize, grid, nu, rho_0, u_0[i]);
    Eigen::ArrayXXd curl = curlZ(state.ux, state.uy, gridsize, dr);

    // ----------- plot the results -----------
    // err = py_figure(pipe, i, std::array<double, 2>{10, 6});
    // if (err != 0) {
    //   print_err("main.cpp", "Failed to execute figure python figure");
    // }
    // err = py_print(pipe, shared_dArray, gridsize, curl);
    // if (err != 0) {
    //   print_err("main.cpp", "Failed to execute print python action");
    // }
    // err = py_pcolor(pipe, shared_dArray, i, gridsize, grid, curl);
    err = py_pcolor(pipe, shared_dArray, semph_ptr, i, gridsize, grid,
                    (state.ux + state.uy).sqrt());
    if (err != 0) {
      print_err("main.cpp", "Failed to execute figure python action");
    }
  }

  err = py_clean(p_sharedmemory_struct, pipe, shared_dArray, semph_ptr);
  if (err != 0) {
    print_err("main.cpp", "Failed to clean up python action");
  }

  return 0;
}
