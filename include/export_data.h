#ifndef EXPORT_DATA_H
#define EXPORT_DATA_H

#include "Eigen/Dense"
#include <map>
#include <string>
#include <variant>

#include "main.h"

typedef std::map<std::string, std::variant<double, std::string>> MetaData;
typedef enum { FORCE, CONFIRM } SaveFlag;

int init_save_dir(const std::string &sim_name, MetaData metadata,
                  SaveFlag save_flag);
int save_state(const std::string &sim_name, State state, Gridsize gridsize,
               Grid grid, double sim_time, SaveFlag save_flag);

#endif // EXPORT_DATA_H
