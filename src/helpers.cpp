#include "helpers.h"
#include <cerrno>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <string.h>

static std::string filename = "helpers.cpp";

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
void print_err(const std::string &filename, const std::string &message) {
  // size_t errmsglen = strerrorlen_s(errno) + 1;
  // size_t errmsglen = 94;
  // char errmsg[errmsglen];
  // strerror_s(errmsg, errmsglen, errno);
  std::cerr << "ERROR <" << filename << ">: [" << strerror(errno)
            << "(error number: " << errno << ")] " << message << std::endl;
}

/**
 * \brief Logs messages to the console or error handler based on the log level.
 *
 * This function logs messages with varying severity levels to the console
 * or an error handler depending on the provided log level. If the log level
 * is greater than the globally defined `LOG_LEVEL`, the message will not be
 * logged.
 *
 * \param[in] filename The name of the file from which the log message
 * originates.
 * \param[in] message The log message to be displayed or handled.
 * \param[in] log_level The severity level of the log message (e.g., error,
 * warning, info, debug).
 *
 * The log levels and their corresponding behavior are:
 * - \c LOG_LEVEL_OFF: No logging.
 * - \c LOG_LEVEL_ERR: Logs errors to the console with the prefix "ERROR" using
 * the `print_err` function.
 * - \c LOG_LEVEL_WRN: Logs warnings to the console with the prefix "WARNING".
 * - \c LOG_LEVEL_INF: Logs informational messages to the console with the
 * prefix "INFO".
 * - \c LOG_LEVEL_DBG: Logs debug messages to the console with the prefix
 * "DEBUG".
 *
 * \see LOG_LEVEL: Global log level.
 * \see Log_level: Different kinds of log levels.
 * \see log: Logging function.
 * \see print_err: Function for printing errors.
 *
 */
void log(const std::string &filename, const std::string &message,
         const Log_level log_level) {
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
