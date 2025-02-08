#include "main.h"
#include "analytical.h"
#include "ex_3.h"
#include "lbm_core.h"

#ifdef GOOGLE_PROF
#include <gperftools/profiler.h>
#endif
#include <string>

static std::string this_filename = "main.cpp";
void initial_condition(State &state, const Gridsize &gridsize, const Grid &grid,
                       [[maybe_unused]] const double nu, const double rho_0,
                       const double u_0, const double p_0);

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
    initCond_TaylorGreen_2D(state, gridsize, grid, rho_0, u_0, p_0);
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

    ex3_main();

#ifdef GOOGLE_PROF
    std::cout << "using goolge profiler" << std::endl;
    ProfilerStop();
#endif
    return 0;
}
