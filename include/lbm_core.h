#ifndef LBM_CORE_H
#define LBM_CORE_H

#include "Eigen/Dense"
#include "ProjectConfig.h"
#include "main.h"

#include <string>

const unsigned int D_ = D;
const unsigned int Q_ = Q;

/**
 * \brief Holds the state variables of the simulation.
 *
 * This structure contains the flow variables: density (rho), velocity
 * components (ux, uy, uz), and pressure (P) for the simulation.
 */
typedef struct {
#if D == 1
    Eigen::ArrayXXd rho; ///< Density array.
    Eigen::ArrayXXd ux;  ///< Velocity component in the x-direction.
    Eigen::ArrayXXd P;   ///< Pressure array.
#elif D == 2
    Eigen::ArrayXXd rho; ///< Density array.
    Eigen::ArrayXXd ux;  ///< Velocity component in the x-direction.
    Eigen::ArrayXXd uy;  ///< Velocity component in the y-direction.
    Eigen::ArrayXXd P;   ///< Pressure array.
#elif D == 3
    Eigen::ArrayXXd rho; ///< Density array.
    Eigen::ArrayXXd ux;  ///< Velocity component in the x-direction.
    Eigen::ArrayXXd uy;  ///< Velocity component in the y-direction.
    Eigen::ArrayXXd uz;  ///< Velocity component in the z-direction.
    Eigen::ArrayXXd P;   ///< Pressure array.
#else
#error "Dimension can't be greater than 3D"
#endif
} State;

/**
 * \brief Specifies the size of the computational grid.
 *
 * This typedef defines an array to hold the dimensions of the computational
 * grid (Nx, Ny, Nz). In case of a dimension lower than 3D the corresponding
 * sizes are simply 1 (for e.g. in 2D: Nz=1).
 */
typedef std::array<std::size_t, 3> Gridsize;

/**
 * \brief Contains the grid vectors for the simulation.
 *
 * This structure holds the 1D grid vectors for the x, y and z coordinates.
 */
typedef struct {
    Eigen::ArrayXd x; ///< Grid vector in the x-direction.
#if D > 1
    Eigen::ArrayXd y; ///< Grid vector in the y-direction.
#endif
#if D > 2
    Eigen::ArrayXd z; ///< Grid vector in the z-direction.
#endif
#if D > 3
#error "Dimension can't be greater than 3D"
#endif
} GridVectors;

/**
 * \brief Contains the grid arrays for the simulation.
 *
 * This structure holds the 2D grid arrays for the X, Y and Z coordinates.
 */
typedef struct {
#if D == 1
    Eigen::ArrayXXd X; ///< 2D grid array in the X-coordinate.
#elif D == 2
    Eigen::ArrayXXd X; ///< 2D grid array in the X-coordinate.
    Eigen::ArrayXXd Y; ///< 2D grid array in the Y-coordinate.
#elif D == 3
    Eigen::ArrayXXd X; ///< 2D grid array in the X-coordinate.
    Eigen::ArrayXXd Y; ///< 2D grid array in the Y-coordinate.
    Eigen::ArrayXXd Z; ///< 2D grid array in the Z-coordinate.
#else
#error "Dimension can't be greater than 3D"
#endif
} Grid;

/**
 * \brief Defines the solver types for the lattice Boltzmann method.
 *
 * This enum lists the different solver types available for the lattice
 * Boltzmann method.
 */
typedef enum {
    LBM,                  ///< Standard LBM solver.
    LBM_EXACT_DIFFERENCE, ///< LBM with exact difference method.
    KBC                   ///< KBC solver.
} SolverType;

int lattice_bolzmann_simulation(
    State &state, const std::string &sim_name, const Gridsize &gridsize,
    const Grid &grid, const double nu, const unsigned int sim_time,
    const SolverType solver,
    std::function<void(State &, const Gridsize &, const Grid &)>
        initial_condition);

#endif // !LBM_CORE_H
