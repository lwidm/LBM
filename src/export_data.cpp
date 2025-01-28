
#include "export_data.h"

#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>
#include <string>

#include "helpers.h"

static std::string this_filename = "export_data.cpp";

int save_eigen_matrix(const Eigen::ArrayXXd &array,
                      const std::string &output_filename, SaveFlag save_flag) {

  if (std::filesystem::exists(output_filename)) {
    LOG_WRN(this_filename,
            "Output state file for current simulation step already exists.");

    if (save_flag != FORCE) {
      char choice;
      std::cout << "The save state file \"" << output_filename
                << "\" already exists. Do you want to overwrite it? (y/n): ";
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

std::string metadata_map_to_json(const MetaData &metadata) {
  std::ostringstream oss;
  oss << "{\n";
  bool first = true;
  for (const auto &[key, value] : metadata) {
    if (!first) {
      oss << ",\n";
    }
    first = false;
    std::visit(
        [&oss](auto &&arg) {
          if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, double>) {
            oss << arg;
          } else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>,
                                              std::string>) {
            oss << "\"" << arg << "\"";
          }
        },
        value);
  }
  oss << "\n}";
  return oss.str();
}

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

int init_save_dir(const std::string &sim_name, MetaData metadata,
                  SaveFlag save_flag) {
  std::filesystem::create_directory(OUTPUT_DIR);
  std::string sim_dir = OUTPUT_DIR + std::string("/") + sim_name;

  if (std::filesystem::exists(sim_dir)) {
    LOG_WRN(this_filename,
            "Output directory for current simulation already exists.");
    if (save_flag == FORCE) {
      std::filesystem::remove_all(sim_dir);
    } else {
      char choice;
      std::cout << "The direcotry \"" << sim_dir
                << "\" already exists. Do you want to remove it? (y/n): ";
      std::cin >> choice;
      if (choice != 'y') {
        std::ostringstream oss;
        oss << "Couldn't save the state of the current simulation: The "
               "direcotry \""
            << sim_dir << "\" already exists.";
        LOG_ERR("this_filename", oss.str());
        return 1;
      }
      std::filesystem::remove_all(sim_dir);
    }
  }
  create_new_save_dir(sim_dir, metadata);
  return 0;
}

int save_state(const std::string &sim_name, State state, Gridsize gridsize,
               Grid grid, double sim_time, SaveFlag save_flag) {

  std::filesystem::create_directory(OUTPUT_DIR);
  std::string sim_dir = OUTPUT_DIR + std::string("/") + sim_name;

  std::ostringstream oss;
  oss << std::scientific << std::setprecision(4) << std::setw(8) << sim_time;
  std::string time_str = oss.str();
  int err;
  err = save_eigen_matrix(state.rho, sim_dir + "/rho_t=" + time_str + ".bin",
                          save_flag);
  if (err != 0) {
    LOG_ERR(this_filename, "Failed to save \"state.rho\".");
  }
  err = save_eigen_matrix(state.ux, sim_dir + "/ux_t=" + time_str + ".bin",
                          save_flag);
  if (err != 0) {
    LOG_ERR(this_filename, "Failed to save \"state.ux\".");
  }
  err = save_eigen_matrix(state.uy, sim_dir + "/uy_t=" + time_str + ".bin",
                          save_flag);
  if (err != 0) {
    LOG_ERR(this_filename, "Failed to save \"state.uy\".");
  }

  return 0;
}
