#include "main.h"
#include "helpers.h"
#include "py_action.h"

#include <array>
#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>
#include <minwindef.h>

void analytical_Poiseuille(State &state, const Gridsize &gridsize,
                           const Grid &grid, const double &nu,
                           const double &rho_0, const double &u_0,
                           const Eigen::ArrayXXd &Fx) {
  state.rho = Eigen::ArrayXXd::Constant(gridsize[1], gridsize[0], rho_0);
  state.uy = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
  state.ux = Fx / (rho_0 * nu) * grid.Y * ((double)(gridsize[1]) - grid.Y);
}

void initial_condition(State &state, const Gridsize &gridsize, const Grid &grid,
                       const double &nu, const double &rho_0,
                       const double &u_0) {
  Eigen::ArrayXXd Fx = Eigen::ArrayXXd::Constant(
      gridsize[1], gridsize[0],
      12.0 * rho_0 * u_0 * nu /
          ((double)(gridsize[1]) * (double)(gridsize[1])));
  analytical_Poiseuille(state, gridsize, grid, nu, rho_0, nu, Fx);
}

int lattice_bolzmann_solver(
    State &state, const Gridsize &gridsize, const Grid &grid,
    const SolverType &solver,
    std::function<void(State &, const Gridsize &, const Grid &)>
        initial_condition) {
  initial_condition(state, gridsize, grid);
  return 0;
}

int main(int argc, char *argv[]) {
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

  // initializing shared array
  SharedMemory_struct sharedmemory_struct;
  std::unique_ptr<double[]> shared_dArray;
  FILE *pipe;
  std::size_t max_size =
      (std::size_t)(22.0 * (double)(D[2]) * 4.1 * (double)(D[2])) + 1;
  std::size_t mmap_max_vector_size = max_size + max_size + max_size;
  int err =
      py_init(sharedmemory_struct, shared_dArray, pipe, mmap_max_vector_size);
  if (err != 0) {
    print_err("main.cpp", "Failed to initialize py_action");
  }

  for (std::size_t i = 0; i < 3; ++i) {

    std::size_t Nx = 22 * D[i];
    std::size_t Ny = (std::size_t)(4.1 * (double)(D[i]));

    const Gridsize gridsize = std::array<std::size_t, 2>{Nx, Ny};

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
    const std::size_t Nz = 1;
    double *curl_plot = curl.data();
    const std::size_t size = Nx * Ny * Nz;
    double *x = gridvectors.x.data();
    double *y = gridvectors.y.data();
    double *z = nullptr;
    double *X = grid.X.data();
    double *Y = grid.Y.data();
    double *Z = nullptr;

    for (std::size_t i = 0; i < size; ++i) {
      shared_dArray[i] = curl_plot[i];
    }

    err = py_print(pipe, shared_dArray, Nx, Ny, Nz);
    if (err != 0) {
      print_err("main.cpp", "Failed to execute print python action");
    }
  }

  err = py_clean(sharedmemory_struct, pipe, shared_dArray);
  if (err != 0) {
    print_err("main.cpp", "Failed to clean up python action");
  }

  return 0;
}
