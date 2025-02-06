#include "analytical.h"
#include "export_data.h"
#include "helpers.h"
#include "lbm_core.h"
#include <array>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <sstream>
#include <string>

const std::string this_filename = "ex_3.cpp";

void ex3_initial_condition(State &state, const Gridsize &gridsize,
                           const Grid &grid, [[maybe_unused]] const double nu,
                           const double rho_0, const double u_0,
                           const double p_0) {
    initCond_TaylorGreen(state, gridsize, grid, rho_0, u_0, p_0);
}

int ex3_LBM_loop(const std::string sim_name, const double dr, const double L,
                 const double t_prime, const double nu, const double rho_0,
                 const double u_0, const double p_0, const SolverType solver) {
    std::size_t Nx = L;
    std::size_t Ny = L;
    std::size_t Nz = 1;
    const unsigned int sim_time = (unsigned int)(t_prime * L / u_0);

    const Gridsize gridsize = std::array<std::size_t, 3>{Nx, Ny, Nz};

    GridVectors gridvectors;
    gridvectors.x = Eigen::ArrayXd::LinSpaced(Nx, dr / 2, Nx * dr);
    gridvectors.y = Eigen::ArrayXd::LinSpaced(Ny, dr / 2, Ny * dr);
    Grid grid = meshgrid(gridsize, gridvectors);

    MetaData metadata;
    std::string description = "A standard LBM simulation of a 2D "
                              "Taylor Green flow. Since the "
                              "analytical solution of a Taylor "
                              "Green flow is know this can be used "
                              "to "
                              "calculate the simulations' error.";
    std::string initial_condition_string = "2D Taylor Green";
    create_metadata(metadata, 1, sim_name, description, gridsize, nu, rho_0,
                    u_0, p_0, solver, initial_condition_string);
    if (init_save_dir(sim_name, metadata, grid, FORCE) != 0) {
        LOG_ERR(this_filename, "Creating new save failed");
    }

    auto initial_condition_preset =
        [&nu, &rho_0, &u_0, &p_0](State &state, const Gridsize &gridsize,
                                  const Grid &grid) {
            ex3_initial_condition(state, gridsize, grid, nu, rho_0, u_0, p_0);
        };

    State state;
    int sim_exit =
        lattice_bolzmann_simulation(state, sim_name, gridsize, grid, nu,
                                    sim_time, solver, initial_condition_preset);
    if (sim_exit != 0) {
        std::ostringstream err_oss;
        err_oss << "Simulation \"" << sim_name << "\" terminated early";
        LOG_WRN(this_filename, err_oss.str());
    }

    Eigen::ArrayXXd curl = curlZ(state.ux, state.uy, gridsize, dr);

    State analytical_state;
    analytical_TaylorGreen(analytical_state, gridsize, grid, nu, rho_0, u_0,
                           p_0, t_prime);
    if (save_state(sim_name, "ana_", analytical_state, (double)sim_time,
                   FORCE) != 0) {
        LOG_ERR(this_filename, "Saving state failed");
        return 1;
    }
    return 0;
}

int ex3_main() {
    // ----------- flow setup -----------
    const double dr = 1;
    // const double dt = 1;
    const double rho_0 = 1;
    const double p_0 = rho_0 / 3;
    const double Re = 100;
    const double fixed_nu = 0.064;
    const double fixed_u_0 = 0.1;
    const double t_prime = 2.5;
    const std::array<std::size_t, 3> L = {64, 96, 128};
    const SolverType solver = LBM;

    // ----------- Run the 3 simulations -----------

    // fixed u_0
    for (std::size_t i = 0; i < 3; ++i) {
        const double nu = fixed_u_0 * L[i] / Re;
        const double u_0 = fixed_u_0;

        std::ostringstream oss;
        oss << "2D_Taylor_Green_fixed_u0_nu=" << std::scientific
            << std::setprecision(4) << std::setw(6) << nu;
        std::string sim_name = oss.str();

        if (ex3_LBM_loop(sim_name, dr, L[i], t_prime, nu, rho_0, u_0, p_0,
                         solver) != 0) {
            LOG_ERR(this_filename, "ex_3 LBM loop failed");
            return 1;
        }
    }

    // fixed nu
    for (std::size_t i = 0; i < 3; ++i) {
        const double u_0 = fixed_nu * Re / L[i];
        const double nu = fixed_nu;

        std::ostringstream oss;
        oss << "2D_Taylor_Green_fixed_nu_u0=" << std::scientific
            << std::setprecision(4) << std::setw(6) << u_0;
        std::string sim_name = oss.str();

        if (ex3_LBM_loop(sim_name, dr, L[i], t_prime, nu, rho_0, u_0, p_0,
                         solver) != 0) {
            LOG_ERR(this_filename, "ex_3 LBM loop failed");
            return 1;
        }
    }

    return 0;
}
