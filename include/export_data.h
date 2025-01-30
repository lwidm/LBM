#ifndef EXPORT_DATA_H
#define EXPORT_DATA_H

#define _USE_MATH_DEFINES

#include "Eigen/Dense"
#include <map>
#include <string>

#include "lbm_core.h"

/**
 * \typedef MetaData
 * \brief Defines a type for metadata storage.
 *
 * This typedef defines a type for storing metadata as a map where each key
 * is a string, and each value is a variant that can be either a double or a
 * string.
 * \todo reference to MetaData creation function
 */
typedef std::map<std::string, std::string> MetaData;

/**
 * \enum SaveFlag
 * \brief Defines the save flag options for saving operations.
 *
 * This enumeration defines the options for the save flag used in saving
 * operations. It allows specifying whether to force overwrite existing files
 * or to confirm before overwriting.
 *
 * - \c FORCE: Force overwrite existing files.
 * - \c CONFIRM: Prompt the user to confirm before overwriting existing files.
 */
typedef enum {
  FORCE,  /**< Force overwrite existing files. */
  CONFIRM /**< Prompt the user to confirm before overwriting existing files. */
} SaveFlag;

int init_save_dir(const std::string &sim_name, const MetaData &metadata,
                  const Grid grid, const SaveFlag save_flag);

int save_state(const std::string &sim_name,
               const std::string &additional_string, const State state,
               const Gridsize gridsize, const Grid grid, const double sim_time,
               const SaveFlag save_flag);

#endif // EXPORT_DATA_H
