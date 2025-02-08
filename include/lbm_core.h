#ifndef LBM_CORE_H
#define LBM_CORE_H

#include "Eigen/Dense"
#include "ProjectConfig.h"
#include "main.h"

#include <string>

/**
 * \brief Holds the state variables of the simulation.
 *
 * This structure contains the flow variables: density (rho), velocity
 * components (ux, uy, uz), and pressure (P) for the simulation.
 */
typedef struct {
    Eigen::ArrayXXd rho; ///< Density array.
    Eigen::ArrayXXd ux;  ///< Velocity component in the x-direction.
#if D > 1
    Eigen::ArrayXXd uy; ///< Velocity component in the y-direction.
#endif
#if D > 2
    Eigen::ArrayXXd uz; ///< Velocity component in the z-direction.
#endif
#if D > 3
#error "Dimension can't be greater than 3D"
#endif
    Eigen::ArrayXXd P; ///< Pressure array.
} State;

/**
 * \brief The descrete velocity set and its corresponding weights.
 *
 * This structure contains the the discrete velocity set separated into its
 * components cxs, cys and czs (depending on the dimension). For e.g. the
 * felocity corresponding to f_3 is (cxs[3],  cys[3], czs[3])^T). The weights
 * used to calculate the equilibruim populations is also saved here
 * */
typedef struct {
    std::array<double, Q> weights;
    std::array<int, Q> cxs;
#if D > 1
    std::array<int, Q> cys;
#endif
#if D > 2
    std::array<int, Q> czs;
#endif
#if D > 3
#error "Dimension can't be greater than 3D"
#endif
} DiscreteVelocities;

/**
 * \brief Specifies the size of the computational grid.
 *
 * This typedef defines an array to hold the dimensions of the computational
 * grid (Nx, Ny, Nz). In case of a dimension lower than 3D the corresponding
 * sizes are simply 1 (for e.g. in 2D: Nz=1).
 */
typedef std::array<Eigen::Index, 3> Gridsize;

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
    Eigen::ArrayXXd X; ///< 2D grid array in the X-coordinate.
#if D > 1
    Eigen::ArrayXXd Y; ///< 2D grid array in the Y-coordinate.
#endif
#if D > 2
    Eigen::ArrayXXd Z; ///< 2D grid array in the Z-coordinate.
#endif
#if D > 3
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
