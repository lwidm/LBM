#include "helpers.h"
#include <cerrno>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <string.h>

static std::string filename = "helpers.cpp";

/**
 * \brief Generates a 2D meshgrid from the provided grid vectors.
 *
 * This function creates 2D meshgrid arrays for the X and Y coordinates based on
 * the provided grid size and grid vectors.
 *
 * \param[in] gridsize The size of the computational grid in Nx and Ny (Nz must
 * be 1).
 * \param[in] gridvectors The grid vectors containing the 1D arrays of x and y
 * coordinates.
 * \return A Grid struct containing the 2D meshgrid arrays for the X and Y
 * coordinates.
 *
 * \details Currently not implemented for 3D Arrays. This function will return
 * an error if Nz is not equal to 1.
 */
Grid meshgrid(const Gridsize &gridsize, const GridVectors &gridvectors) {
  const std::size_t Nx = gridsize[0];
  const std::size_t Ny = gridsize[1];
  if (gridsize[2] != 1) {
    LOG_ERR(filename, "meshgrid function only implemented for 2D grids");
    return Grid();
  }
  const Eigen::ArrayXd x = gridvectors.x;
  const Eigen::ArrayXd y = gridvectors.y;

  Eigen::ArrayXXd X(Ny, Nx);
  Eigen::ArrayXXd Y(Ny, Nx);

  // TODO : Use the Eigen row assignements for speed
  for (std::size_t i = 0; i < Nx; ++i) {
    for (std::size_t j = 0; j < Ny; ++j) {
      X(j, i) = x(i);
      Y(j, i) = y(j);
    }
  }

  Grid grid = {X, Y};
  return grid;
}

/**
 * \brief Calculates the z-component of the curl for a 2D velocity field.
 *
 * This function computes the z-component of the curl (vorticity) for a given 2D
 * velocity field represented by the ux and uy components. It uses central
 * differences for the numerical derivative calculation.
 *
 * \param[in] ux The velocity component in the x-direction.
 * \param[in] uy The velocity component in the y-direction.
 * \param[in] gridsize The size of the computational grid in Nx and Ny (Nz must
 * be 1).
 * \param[in] dr The grid spacing (assumed to be equal in both directions).
 * \return An Eigen::ArrayXXd containing the z-component of the curl for the 2D
 * velocity field.
 *
 * \details This function will return an error if Nz is not equal to 1.
 * \details Currently, this function does not handle boundary conditions
 * properly.
 */
Eigen::ArrayXXd curlZ(const Eigen::ArrayXXd &ux, const Eigen::ArrayXXd &uy,
                      const Gridsize &gridsize, const double dr) {
  const std::size_t Nx = gridsize[0];
  const std::size_t Ny = gridsize[1];
  if (gridsize[2] != 1) {
    LOG_ERR(filename, "curlZ does only supports 2D grids");
    return Eigen::ArrayXXd();
  }
  const double dx = dr;
  const double dy = dr;
  // BUG : Do this propery with the boundery condition (Matrix operation)
  Eigen::ArrayXXd curl = Eigen::ArrayXXd::Constant(Ny, Nx, 0);
  for (std::size_t i = 1; i < Nx - 1; ++i) {
    for (std::size_t j = 1; j < Ny - 1; ++j) {
      double duydx = (uy(j, i + 1) - uy(j, i - 1)) / (2 * dx);
      double duxdy = (ux(j + 1, i) - ux(j - 1, i)) / (2 * dy);
      curl(j, i) = duydx - duxdy;
    }
  }
  return curl;
}

/**
 * \brief Prints an error message to the standard output.
 *
 * This function prints an error message to the standard output, including the
 * file name, a description of the error based on the current value of `errno`,
 * and a custom message.
 *
 * \param[in] filename The name of the file where the error occurred.
 * \param[in] message A custom message describing the context or details of the
 * error.
 */
void print_err(std::string filename, std::string message) {
  // size_t errmsglen = strerrorlen_s(errno) + 1;
  // size_t errmsglen = 94;
  // char errmsg[errmsglen];
  // strerror_s(errmsg, errmsglen, errno);
  std::cerr << "ERROR <" << filename << ">: [" << strerror(errno)
            << "(error number: " << errno << ")] " << message << std::endl;
}

void log(std::string filename, std::string message, Log_level log_level) {
  if (log_level > LOG_LEVEL) {
    return;
  }
  switch (log_level) {
  case LOG_LEVEL_OFF:
    break;
  case LOG_LEVEL_ERR:
    print_err(filename, message);
    break;
  case LOG_LEVEL_WRN:
    std::cout << "WARNING <" << filename << ">: " << message << std::endl;
    break;
  case LOG_LEVEL_INF:
    std::cout << "INFO <" << filename << ">: " << message << std::endl;
    break;
  case LOG_LEVEL_DBG:
    std::cout << "DEBUG <" << filename << ">: " << message << std::endl;
    break;
  }
}
