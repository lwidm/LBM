
#include "export_data.h"
#include "ProjectConfig.h"

#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>
#include <string>

#include "helpers.h"

static std::string this_filename = "export_data.cpp";

/**
 * \brief Saves a 2D Eigen matrix to a binary file.
 *
 * This function saves an Eigen matrix to a specified binary file. If the file
 * already exists and the save flag is not set to FORCE, the user is prompted
 * to confirm overwriting the file.
 *
 * \param[in] array The Eigen matrix to be saved.
 * \param[in] output_filename The name of the output file.
 * \param[in] save_flag The flag indicating whether to force overwrite.
 *
 * \return 0 if the matrix is saved successfully, 1 if an error occurs.
 *
 * **Warning**: Eigen uses column-major ordering by default so this function
 * saves using column-major ordering. Many programs (including NumPy in
 * Python) use row-major ordering by default.
 *
 * \see load_eigen_matrix: Python function to load an Eigen matrix from a binary
 * file
 */
int save_eigen_matrix(const Eigen::ArrayXXd &array,
                      const std::string &output_filename,
                      const SaveFlag save_flag) {

  if (std::filesystem::exists(output_filename)) {
    LOG_WRN(this_filename,
            "Output state file for current simulation step already exists.");

    if (save_flag != FORCE) {
      char choice;
      std::cout << "The save state file \"" << output_filename
                << "\" already exists.\nDo you want to overwrite it? (y/n): ";
      std::cin >> choice;

      if (choice != 'y') {
        std::ostringstream oss;
        oss << "Couldn't save the state of the current simulation step: The "
               "file \""
            << output_filename << "\" already exists.";
        LOG_ERR("this_filename", oss.str());
        return 1;
      }
    }
  }

  std::ofstream outFile(output_filename, std::ios::binary);
  if (outFile.is_open()) {
    unsigned int rows = array.rows();
    unsigned int cols = array.rows();
    outFile.write((char *)&rows, sizeof(unsigned int));
    outFile.write((char *)&cols, sizeof(unsigned int));
    outFile.write((char *)array.data(), rows * cols * sizeof(double));
    outFile.close();
  } else {
    LOG_ERR(this_filename,
            "Unable to open file in function: \"save_eigen_matrix(...)\"");
    return 0;
  }
  return 0;
}

/**
 * \brief Converts metadata to a JSON string.
 *
 * This function converts a metadata map to a JSON string format.
 *
 * \param[in] metadata The metadata map to be converted.
 *
 * \return A JSON string representation of the metadata.
 */
std::string metadata_map_to_json(const MetaData &metadata) {
  std::ostringstream oss;
  oss << "{\n";
  bool first = true;
  for (const auto &[key, value] : metadata) {
    if (!first) {
      oss << ",\n";
    } else {
      first = false;
    }
    oss << "\"" << key << "\": " << value;
  }
  oss << "\n}";
  return oss.str();
}

/**
 * \brief Creates a new save directory and saves the simulations metadata to it.
 *
 * This function creates a new save directory for the simulation and saves
 * the metadata as a JSON file in the directory.
 *
 * \param[in] sim_dir The name of the simulation directory.
 * \param[in] metadata The metadata to be saved in the directory.
 *
 * \return 0 if the directory and metadata file are created successfully, 1 if
 * an error occurs.
 */
int create_new_save_dir(const std::string &sim_dir, const MetaData &metadata) {
  std::filesystem::create_directory(sim_dir);
  std::string metadata_json = metadata_map_to_json(metadata);
  std::string metadata_filename = sim_dir + "/metadata.json";
  std::ofstream metadata_file(metadata_filename);
  if (metadata_file.is_open()) {
    metadata_file << metadata_json;
    metadata_file.close();
  } else {
    LOG_ERR(this_filename, "Unable to open metadata json file");
    return 1;
  }
  return 0;
}

/**
 * \brief Initializes the save directory for the simulation.
 *
 * This function initializes the save directory for a simulation. If the
 * directory already exists and the save flag is not set to FORCE, the user is
 * prompted to confirm removal of the existing directory.
 *
 * \param[in] sim_name The name of the simulation.
 * \param[in] metadata The metadata for the simulation.
 * \param[in] grid A Grid struct containing the 2D meshgrid arrays for the X (,Y
 * and Z) coordinates.
 * \param[in] save_flag The flag indicating whether to force remove existing
 * directory.
 *
 * \return 0 if the directory is initialized successfully, 1 if an error occurs.
 */
int init_save_dir(const std::string &sim_name, const MetaData &metadata,
                  const Grid grid, const SaveFlag save_flag) {
  std::filesystem::create_directory(OUTPUT_DIR);
  std::string sim_dir = OUTPUT_DIR + std::string("/") + sim_name;

  if (std::filesystem::exists(sim_dir)) {
    LOG_WRN(this_filename,
            "Output directory for current simulation already exists.");
    if (save_flag != FORCE) {
      char choice;
      std::cout << "The direcotry \"" << sim_dir
                << "\" already exists.\nDo you want to remove it? (y/n): ";
      std::cin >> choice;
      if (choice != 'y') {
        std::ostringstream oss;
        oss << "Couldn't save the state of the current simulation: The "
               "direcotry \""
            << sim_dir << "\" already exists.";
        LOG_ERR("this_filename", oss.str());
        return 1;
      }
    }
    std::filesystem::remove_all(sim_dir);
  }
  create_new_save_dir(sim_dir, metadata);
  save_eigen_matrix(grid.X, sim_dir + "/X.bin", save_flag);
  save_eigen_matrix(grid.Y, sim_dir + "/Y.bin", save_flag);
  return 0;
}

/**
 * \brief Saves the state of the simulation.
 *
 * \param[in] sim_name The name of the simulation.
 * \param[in] additional_string A string that gets added to the start of the
 * filenames (often one leaves this empty "")
 * \param[in] state The state of the simulation to be saved.
 * \param[in] sim_time The simulation time at which the state is saved.
 * \param[in] save_flag The flag indicating whether to force overwrite existing
 * files.
 *
 * \return 0 if the state is saved successfully, 1 if an error occurs.
 *
 * This function saves the state of the simulation at a specific time step.
 * The state includes matrices for density (rho) and velocity components (ux,
 * uy). It saves files in the format
 * \<additional_string\>\<state_variable\>_t=X.XXXXe+XX.bin
 *
 * Example:
 * - additional_string = "analytical_"
 * - sim_time = 200
 * => analytical_ux_t=2.0000e+02.bin, analytical_uy_t=2.0000e+02.bin, etc.
 */
int save_state(const std::string &sim_name,
               const std::string &additional_string, const State state,
               const double sim_time, const SaveFlag save_flag) {

  std::filesystem::create_directory(OUTPUT_DIR);
  std::string sim_dir = OUTPUT_DIR + std::string("/") + sim_name;

  std::ostringstream oss;
  oss << std::scientific << std::setprecision(6) << std::setw(10) << sim_time;
  std::string time_str = oss.str();
  int err;
  err = save_eigen_matrix(state.rho,
                          sim_dir + "/" + additional_string +
                              "rho_t=" + time_str + ".bin",
                          save_flag);
  if (err != 0) {
    LOG_ERR(this_filename, "Failed to save \"state.rho\".");
  }
  err = save_eigen_matrix(
      state.ux, sim_dir + "/" + additional_string + "ux_t=" + time_str + ".bin",
      save_flag);
  if (err != 0) {
    LOG_ERR(this_filename, "Failed to save \"state.ux\".");
  }
  err = save_eigen_matrix(
      state.uy, sim_dir + "/" + additional_string + "uy_t=" + time_str + ".bin",
      save_flag);
  if (err != 0) {
    LOG_ERR(this_filename, "Failed to save \"state.uy\".");
  }

  return 0;
}

/**
 * \brief Creates a metadata object for a simulation and stores it in a map.
 *
 * This function generates metadata for a simulation and stores it in a
 * `MetaData` map. The metadata includes essential information such as the
 * simulation ID, name, description and other parameters
 *
 * The function uses `std::ostringstream` to convert numerical values (e.g.,
 * `double`, `int`) to strings with appropriate precision and formatting.
 *
 * \param[in,out] metadata The map where the metadata will be stored.
 * \param[in] id The unique identifier for the simulation.
 * \param[in] sim_name The name of the simulation.
 * \param[in] description A brief description of the simulation.
 * \param[in] gridsize
 * \param[in] nu
 * \param[in] rho_0
 * \param[in] u_0
 * \param[in] p_0
 * \param[in] solver
 * \param[in] initial_condition
 *
 * \see MetaData: The map type to store metadata
 */
void create_metadata(MetaData &metadata, const unsigned int id,
                     const std::string &sim_name,
                     const std::string &description, const Gridsize &gridsize,
                     const double nu, const double rho_0, const double u_0,
                     const double p_0, const SolverType solver,
                     const std::string initial_condition) {

  std::ostringstream oss_int;

  std::ostringstream oss_double;
  oss_double << std::setprecision(8) << std::defaultfloat;

  oss_int << id;
  metadata["id"] = oss_int.str();
  oss_int.clear();
  oss_int.str("");

  metadata["sim_name"] = "\"" + sim_name + "\"";

  metadata["description"] = "\"" + description + "\"";

  oss_int << gridsize[0];
  metadata["Nx"] = oss_int.str();
  oss_int.clear();
  oss_int.str("");

  oss_int << gridsize[1];
  metadata["Ny"] = oss_int.str();
  oss_int.clear();
  oss_int.str("");

  oss_int << gridsize[2];
  metadata["Nz"] = oss_int.str();
  oss_int.clear();
  oss_int.str("");

  oss_double << nu;
  metadata["nu"] = oss_double.str();
  oss_double.clear();
  oss_double.str("");

  oss_double << rho_0;
  metadata["rho_0"] = oss_double.str();
  oss_double.clear();
  oss_double.str("");

  oss_double << u_0;
  metadata["u_0"] = oss_double.str();
  oss_double.clear();
  oss_double.str("");

  oss_double << p_0;
  metadata["p_0"] = oss_double.str();
  oss_double.clear();
  oss_double.str("");

  std::string solver_string;
  switch (solver) {
  case LBM:
    solver_string = "LBM";
    break;
  case LBM_EXACT_DIFFERENCE:
    solver_string = "LBM exact difference";
    break;
  case KBC:
    solver_string = "KBC";
    break;
  }
  metadata["solver"] = solver_string;

  metadata["initial condition"] = initial_condition;
}
