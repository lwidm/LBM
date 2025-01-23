#include "helpers.h"
#include <cerrno>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <string.h>

Grid meshgrid(const Gridsize &gridsize, const GridVectors &gridvectors) {
  const std::size_t Nx = gridsize[0];
  const std::size_t Ny = gridsize[1];
  const Eigen::ArrayXd x = gridvectors.x;
  const Eigen::ArrayXd y = gridvectors.y;

  Eigen::ArrayXXd X(Ny, Nx);
  Eigen::ArrayXXd Y(Ny, Nx);

  for (std::size_t i = 0; i < Nx; ++i) {
    for (std::size_t j = 0; j < Ny; ++j) {
      Y(j, i) = y(j);
      X(j, i) = x(i);
    }
  }

  Grid grid = {X, Y};
  return grid;
}

Eigen::ArrayXXd curlZ(const Eigen::ArrayXXd &ux, const Eigen::ArrayXXd &uy,
                      const Gridsize &gridsize, const double dr) {
  const std::size_t Nx = gridsize[0];
  const std::size_t Ny = gridsize[1];
  const double dx = dr;
  const double dy = dr;
  Eigen::ArrayXXd curl(Ny, Nx);
  for (std::size_t i = 1; i < Nx - 1; ++i) {
    for (std::size_t j = 1; j < Ny - 1; ++j) {
      double duydx = (uy(j, i + 1) - uy(j, i - 1)) / (2 * dx);
      double duxdy = (ux(j + 1, i) - ux(j - 1, i)) / (2 * dy);
      curl(j, i) = duydx - duxdy;
    }
  }
  return curl;
}

void print_err(std::string filename, std::string message) {
  // size_t errmsglen = strerrorlen_s(errno) + 1;
  // size_t errmsglen = 94;
  // char errmsg[errmsglen];
  // strerror_s(errmsg, errmsglen, errno);
  std::cout << "ERROR <" << filename << ">: [" << strerror(errno)
            << "(error number: " << errno << ")] " << message;
}
