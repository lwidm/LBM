#include "lbm_core.h"
#include "export_data.h"
#include "helpers.h"

#include <array>
#include <iostream>

const std::string this_filename = "lbm_core.cpp";
void streaming_step(std::array<Eigen::ArrayXXd, Q> &f, const Gridsize &gridsize,
                    const DiscreteVelocities &discreteVelocities);
void compute_macroscopic_variables(State &state, const Gridsize &gridsize,
                                   const std::array<Eigen::ArrayXXd, Q> &f,
                                   const DiscreteVelocities &discreteVelocities,
                                   const double cs);

/**
 * \brief The streaming step of the Lattice Baltzmann method
 *
 * This function applies the streaming to the populations f. To determine in
 * which direction every population f[i] should be streamed the vecotors csx and
 * csy are used. f cxs and cys are array of length Q where Q is the number of
 * populations in the simulation (e.g. D2Q9 -> Q = 9, D1Q3 -> Q = 3). The
 * cirshift function is used for this computation.
 *
 * \param[in,out] f The populations to be shifted
 * \param[in] gridsize The size of the grid in the format [Nx, Ny, Nz] (Nz = 1)
 * \param[in] discreteVelocities The struct containing the discrete velocity
 * components of the Lattice [cxs, cys, czs].
 *
 * \see circshift: Cirshift operation on Eigen Arrays
 */
void streaming_step(std::array<Eigen::ArrayXXd, Q> &f, const Gridsize &gridsize,
                    const DiscreteVelocities &discreteVelocities) {
    for (unsigned int i = 1; i < Q; ++i) {
#if D == 1
#error "1D roll function not implemented yet"
#elif D == 2
        roll2D(f[i], gridsize, discreteVelocities.cys[i],
               discreteVelocities.cxs[i]);
#elif D == 3
#error "3D roll function not implemented yet"
#else
#error "Dimension can't be greater than 3D"
#endif
    }
}

/**
 * \brief Compute macroscopic variables (e.g. rho, ux, uy, uz, P) from the
 * populations f
 *
 * This function computes the macroscopic variables rho, ux, uy and uz from the
 * populations f
 *
 * \param[in,out] state The state of the simulations which is a struct that
 * hold the arrays for rho, ux, uy, uz and P (i.e. the state)
 * \param[in] f The current population distribution of the simulation
 * \param[in] gridsize The size of the grid in the format [Nx, Ny, Nz]
 * \param[in] discreteVelocities The struct containing the discrete velocity
 * components of the Lattice [cxs, cys, czs].
 * \param[in] cs The lattice speed of sound
 */
void compute_macroscopic_variables(State &state, const Gridsize &gridsize,
                                   const std::array<Eigen::ArrayXXd, Q> &f,
                                   const DiscreteVelocities &discreteVelocities,
                                   const double cs) {
    state.rho = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
    state.P = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
#if D == 1
    state.ux = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
    for (unsigned int i = 0; i < Q; ++i) {
        state.rho += f[i];
        state.ux += cxs[i] * f[i];
        state.P += cxs[i] * cxs[i] * f[i]
    }
#elif D == 2
    state.ux = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
    state.uy = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
    for (unsigned int i = 0; i < Q; ++i) {
        state.rho += f[i];
        state.ux += discreteVelocities.cxs[i] * f[i];
        state.uy += discreteVelocities.cys[i] * f[i];
    }
    state.ux = state.ux / state.rho;
    state.uy = state.uy / state.rho;
    // BUG : Not sure about pressure term
    state.P = cs * cs + (state.ux * state.ux + state.uy * state.uy);
}
#elif D == 3
    state.ux = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
    state.uy = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
    state.uz = Eigen::ArrayXXd::Zero(gridsize[1], gridsize[0]);
    for (unsigned int i = 0; i < Q; ++i) {
        state.rho += f[i];
        state.ux += cxs[i] * f[i];
        state.uy += cys[i] * f[i];
        state.uz += czs[i] * f[i];
    }
    state.ux = state.ux / state.rho;
    state.uy = state.uy / state.rho;
    state.uz = state.uz / state.rho;
    // BUG : Not sure about pressure term
    state.P = cs * cs +
              (state.ux * state.ux + state.uy * state.uy + state.uz * state.uz);
}
#else
#error "Dimension can't be greater than 3D"
#endif

    /**
     * \brief Main simulation function for the lattice Boltzmann method.
     *
     * This function runs the lattice Boltzmann simulation for the specified
     * simulation time. It modifies the state variables rho, ux, (uy, uz) and P
     * according to the chosen solver type. The function initializes the state
     * using the provided initial condition function and then proceeds with the
     * simulation.
     *
     * \param[in,out] state The state object that holds the flow variables rho,
     * ux, (uy, uz) and P.
     * \param[in] sim_name The name of the current simulation (used for
     * exporting data)
     * \param[in] gridsize The size of the computational grid in Nx, Ny, Nz (Nz
     * = 1 for 2D simulations).
     * \param[in] grid The grid object containing coordinate information
     * (1D/2D/3D arrays of X, X,Y, X,Y,Z depending on the simulation dimension).
     * \param[in] nu The viscosity of the fluid
     * \param[in] sim_time The total simulation time for which the solver should
     * run.
     * \param[in] solver An enum representing the type of lattice Boltzmann
     * solver to be used (e.g., LBM, LBM_EXACT_DIFFERENCE, or KBC). solver to be
     * used (e.g., "LBM", "KBC", or "LBM_exact_difference").
     * \param[in] initial_condition A function that sets the initial conditions
     * for the simulation.
     * \return true if the entire simulation was run for the specified sim_time,
     * false if it terminated early.
     *
     * \details The function returns false if the simulation terminates early.
     * Early termination are implementation specific and depend on the users
     * implementation. For e.g. one may want the function to terminate early if
     * the simulation becomes unstable or convergences.
     */
    int lattice_bolzmann_simulation(
        State & state, const std::string &sim_name, const Gridsize &gridsize,
        const Grid &grid, const double nu, const unsigned int sim_time,
        [[maybe_unused]] const SolverType solver,
        std::function<void(State &, const Gridsize &, const Grid &)>
            initial_condition) {

        // start log
        std::cout << "Starting simulation: " << sim_name << std::endl;

        // output definitions
        const unsigned int output_freq = 20;
        const unsigned int log_freq = 200;
        const bool save_output = true;

        // D2Q9 simulation
        const unsigned int dr =
            1; // BUG : Not sure if this should depend on grid
        const unsigned int dt = 1; // BUG : Not sure if this should be an input
        // const double c = (double)dr / (double)dt;
        const int c = 1;
        const double cs = c / std::sqrt(3.);
        const double beta = static_cast<double>(dt) /
                            (2 * nu / (cs * cs) + static_cast<double>(dt));

        DiscreteVelocities discreteVelocities;
#if D == 1
        discreteVelocities.cxs = {-c, 0, c};
        discreteVelocities.weigths = {1. / 6., 2. / 3., 1. / 6.};
#elif D == 2
    discreteVelocities.cxs = std::array<int, Q>{0, 0, c, c, c, 0, -c, -c, -c};
    discreteVelocities.cys = std::array<int, Q>{0, c, c, 0, -c, -c, -c, 0, c};
    discreteVelocities.weights =
        std::array<double, Q>{4. / 9, 1. / 9,  1. / 36, 1. / 9, 1. / 36,
                              1. / 9, 1. / 36, 1. / 9,  1. / 36};
#elif D == 3
#error "3D lattice_bolzmann_simulation function not implemented yet"
#else
#error "Dimension can't be greater than 3D"
#endif

        auto feq_calc = [&discreteVelocities,
                         &cs](std::array<Eigen::ArrayXXd, Q> &feq,
                              const State &state_) {
#if D == 1
            feq[1] = state_.rho * (1. - (cs * cs + state_.ux * state_.ux));
            feq[2] =
                1. / 2. * state_.rho * (-u + (cs * cs + state_.ux * state_.ux));
            feq[3] =
                1. / 2. * state_.rho * (u + (cs * cs + state_.ux * state_.ux));
#elif D == 2
        for (unsigned int i = 0; i < Q; ++i) {
            std::array<double, Q> weights = discreteVelocities.weights;
            std::array<int, Q> cxs = discreteVelocities.cxs;
            std::array<int, Q> cys = discreteVelocities.cys;
            feq[i] =
                state_.rho * weights[i] *
                (1.0 + 3. * (cxs[i] * state_.ux + cys[i] * state_.uy) +
                 9. * (cxs[i] * state_.ux + cys[i] * state_.uy).square() / 2 -
                 3 * (state_.ux.square() + state_.uy.square()) / 2);
        }
#elif D == 3
#error "3D lattice_bolzmann_simulation function not implemented yet"
#else
#error "Dimension can't be greater than 3D"
#endif
        };

        // Initial condition
        initial_condition(state, gridsize, grid);
        std::array<Eigen::ArrayXXd, Q> f;
        std::array<Eigen::ArrayXXd, Q> feq;
        std::array<Eigen::ArrayXXd, Q> fmirr;
        feq_calc(feq, state);
        f = feq;

        if (save_output) {
            if (save_state(sim_name, "num_", state, 0., CONFIRM) != 0) {
                LOG_ERR(this_filename, "Saving state failed");
                return 1;
            }
        }

        for (unsigned int t = 1; t <= sim_time; ++t) {
            // Collision step
            feq_calc(feq, state);
            for (unsigned int i = 0; i < Q; ++i) {
                fmirr[i] = 2 * feq[i] - f[i];
                f[i] = (1 - beta) * f[i] + beta * fmirr[i];
            }
            // Streaming step
            streaming_step(f, gridsize, discreteVelocities);
            // Boundary Condition
            // Compute macroscopic variables
            compute_macroscopic_variables(state, gridsize, f,
                                          discreteVelocities, cs);

            if (save_output && t % output_freq == 0) {
                if (save_state(sim_name, "num_", state, static_cast<double>(t),
                               CONFIRM) != 0) {
                    LOG_ERR(this_filename, "Saving state failed");
                    return 1;
                }
            }
            if (t % log_freq == 0) {
                std::cout << "\r time: " << t << "/" << sim_time << std::flush;
            }
        }

        // end log
        std::cout << "\n" << std::endl;

        return 0;
    }
