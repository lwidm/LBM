import numpy as np
import matplotlib.pyplot as plt
from enum import Enum
import re
import glob
import time

# from numpy.typing import NDArray
from typing import List

arrayXXXd = np.ndarray[tuple[int, int, int], np.dtype[np.float64]]
arrayXXd = np.ndarray[tuple[int, int], np.dtype[np.float64]]
arrayXd = np.ndarray[tuple[int], np.dtype[np.float64]]


## Defines the various log levels for logging messages.
#
# This enumeration defines different levels of logging messages that can be used to categorize the importance or type of message being logged.
#
# \param LOG_LEVEL_OFF (int): No logging.
# \param LOG_LEVEL_ERR (int): Logs errors to the console with the prefix "ERROR".
# \param LOG_LEVEL_WRN (int): Logs warnings to the console with the prefix "WARNING".
# \param LOG_LEVEL_INF (int): Logs informational messages to the console with the prefix "INFO".
# \param LOG_LEVEL_DBG (int): Logs debug messages to the console with the prefix "DEBUG".
# fmt: off
class Log_level(Enum):
    LOG_LEVEL_OFF = 0  ##< No logging.
    LOG_LEVEL_ERR = 1  ##< Logs errors to the console with the prefix "ERROR".
    LOG_LEVEL_WRN = 2  ##< Logs warnings to the console with the prefix "WARNING".
    LOG_LEVEL_INF = 3  ##< Logs informational messages to the console with the prefix "INFO".
    LOG_LEVEL_DBG = 4  ##< Logs debug messages to the console with the prefix "DEBUG".
# fmt: on


# %% %%%%%%%%%%%%%%%%%%%% Definitions ####################

## Defines the global log level for the logging system.
#
# This variable defines the global log level for the logging system. Messages with a log level higher than this value will not be logged. The log level can be adjusted to control the verbosity of logging output.
#
# \see log (function): Logging function.
# \see Log_level (Enum): Different kinds of log levels.
LOG_LEVEL: Log_level = Log_level.LOG_LEVEL_DBG


this_filename: str = "post_processing.py"


# %% %%%%%%%%%%%%%%%%%%%% Logging functions ####################


## \brief Logs messages to the console based on the log level.
# \param filename (str): The name of the file from which the log message originates.
# \param message (str): The log message to be displayed.
# \param log_level (Log_level): The severity level of the log message (e.g., error, warning, info, debug).
#
# \return None
#
# Log levels and their corresponding behavior:
#   - \c LOG_LEVEL_OFF: No logging.
#   - \c LOG_LEVEL_ERR: Logs errors to the console with the prefix "ERROR".
#   - \c LOG_LEVEL_WRN: Logs warnings to the console with the prefix "WARNING".
#   - \c LOG_LEVEL_INF: Logs informational messages to the console with the prefix "INFO".
#   - \c LOG_LEVEL_DBG: Logs debug messages to the console with the prefix "DEBUG".
#
# \see LOG_LEVEL (Log_level): Global log level.
# \see Log_level (Enum): Different kinds of log levels.
def log(filename: str, message: str, log_level: Log_level) -> None:
    if log_level.value > LOG_LEVEL.value:
        return
    match log_level:
        case Log_level.LOG_LEVEL_ERR:
            print(f"ERROR <{filename}>: {message}")
        case Log_level.LOG_LEVEL_WRN:
            print(f"WARNING <{filename}>: {message}")
        case Log_level.LOG_LEVEL_INF:
            print(f"INFO <{filename}>: {message}")
        case Log_level.LOG_LEVEL_DBG:
            print(f"DEBUG <{filename}>: {message}")
        case _:
            print(f"ERROR <{filename}>: Invalid log level.")
    return


## \brief Logs a error message.
#
# This function logs a error message using the `log` function with the error log level.
#
# \param filename (str): The name of the file from which the log error originates.
# \param message (str): The warning message to be logged.
#
# \see LOG_LEVEL (Log_level): Globally set log-level.
# \see Log_level (Enum): Different kinds of log levels.
# \see log (function): Logging function.
def LOG_ERR(filename: str, message: str) -> None:
    log(filename, message, Log_level.LOG_LEVEL_ERR)
    return


## \brief Logs a warning message.
#
# This function logs a warning message using the `log` function with the warning log level.
#
# \param filename (str): The name of the file from which the log message originates.
# \param message (str): The warning message to be logged.
#
# \see LOG_LEVEL (Log_level): Globally set log-level.
# \see Log_level (Enum): Different kinds of log levels.
# \see log (function): Logging function.
def LOG_WRN(filename: str, message: str) -> None:
    log(filename, message, Log_level.LOG_LEVEL_WRN)
    return


## \brief Logs a info message.
#
# This function logs a info message using the `log` function with the info log level.
#
# \param filename (str): The name of the file from which the log message originates.
# \param message (str): The info message to be logged.
#
# \see LOG_LEVEL (Log_level): Globally set log-level.
# \see Log_level (Enum): Different kinds of log levels.
# \see log (function): Logging function.
def LOG_INF(filename: str, message: str) -> None:
    log(filename, message, Log_level.LOG_LEVEL_INF)
    return


## \brief Logs a debug message.
#
# This function logs a debug message using the `log` function with the debug log level.
#
# \param filename (str): The name of the file from which the log message originates.
# \param message (str): The debug message to be logged.
#
# \see LOG_LEVEL (Log_level): Globally set log-level.
# \see Log_level (Enum): Different kinds of log levels.
# \see log (function): Logging function.
def LOG_DBG(filename: str, message: str) -> None:
    log(filename, message, Log_level.LOG_LEVEL_DBG)
    return


# %% %%%%%%%%%%%%%%%%%%%% Loading matrices from files functions ####################


## \brief Loads an Eigen vector from a binary file.
#
# This function loads an Eigen vector from a specified binary file.
#
# \param filename (str): The name of the file from which the Eigen vector is to be loaded.
#
# \returns arrayXd: The loaded Eigen vector.
def load_eigen_vector(filename: str) -> arrayXd:
    vector: arrayXd
    with open(filename, "rb") as f:
        vector = np.fromfile(f, dtype=np.float64)
    return vector


##
# \brief Loads an Eigen matrix from a binary file.
# \param filename (str): The name of the file from which the Eigen matrix is to be loaded.
# \param order (str, optional): The order in which the matrix should be reshaped. "C" for C-style (row-major), "F" for Fortran-style (column-major). Defaults to "F".
#
# \returns (arrayXXd): The loaded Eigen matrix.
# \exception ValueError: If the 'order' argument is not "C" or "F".
#
# **Warning**: Eigen (the library used for arrays in the c++ program) uses column-major ordering by default.
#
# \see save_eigen_matrix: c++ Function to save Eigen matrix (in column-major form)
def load_eigen_matrix(filename: str, order: str = "F") -> arrayXXd:
    matrix: arrayXXd
    with open(filename, "rb") as f:
        rows: int = np.fromfile(f, dtype=np.uint32, count=1)[0]
        cols: int = np.fromfile(f, dtype=np.uint32, count=1)[0]
        data: arrayXd = np.fromfile(f, dtype=np.float64)
        matrix = data.reshape((rows, cols))
        match order:
            case "C":
                pass
            case "F":
                matrix = matrix.reshape((cols, rows)).T
            case _:
                log(
                    filename,
                    'Invalid "order" value in function "load_eigen_array(filename, order)".',
                    Log_level.LOG_LEVEL_ERR,
                )

    return matrix


# %% %%%%%%%%%%%%%%%%%%%% main block ####################


##
# \brief get a floating point time value from the data filename
#
# Extracts a floating point time value from the time values contained in the datas filename. E.g. "ux_t=03.300e+02.bin" would return "330"
#
# \param filename (str): The string of the filename to be processed
def extract_time_from_filename(filename: str) -> float:
    match = re.search(r"t=(\d+\.\d+e[+-]\d+)", filename)
    time: float
    if match:
        time = float(match.group(1))
    else:
        time = -1
        LOG_ERR(this_filename, f'No time string found in filename: "{filename}".')
    return time


##
# \brief Entry point of the `post_precessing.py` file
def main() -> int:
    directory: str = "../../output/taylor green"

    ux_files: List[str] = sorted(
        glob.glob(directory + "/num_ux_t=*.bin"), key=extract_time_from_filename
    )
    uy_files: List[str] = sorted(
        glob.glob(directory + "/num_uy_t=*.bin"), key=extract_time_from_filename
    )

    if len(ux_files) != len(uy_files):
        LOG_ERR(
            this_filename,
            f"the number of ux_files ({len(ux_files)}) and uy_files ({len(uy_files)}) aren't equal.",
        )
        return 1
    if len(ux_files) == 0:
        LOG_ERR(this_filename, "No ux_files and uy_files found.")
        return 1

    max_time: float = extract_time_from_filename(ux_files[-1])

    ux: arrayXXd
    uy: arrayXXd
    max_mag: float = 0
    plt.ion()
    fig = plt.figure()
    for i in np.arange(len(ux_files)):
        plt.clf()
        if not plt.fignum_exists(fig.number):
            LOG_INF(this_filename, "plot closed before all states were plotted")
            return 0
        ux = load_eigen_matrix(ux_files[i])
        uy = load_eigen_matrix(uy_files[i])
        t: float = extract_time_from_filename(ux_files[i])
        mag = np.sqrt(ux**2, uy**2)
        max_mag = np.max([max_mag, np.max(mag)])
        plt.pcolor(mag)
        plt.title(f"time = {t}/{max_time}")
        plt.clim(0, max_mag)
        plt.colorbar()
        plt.draw()
        plt.pause(0.05)
    ux = load_eigen_matrix(ux_files[-1])
    uy = load_eigen_matrix(uy_files[-1])

    ux_file_ana: str = glob.glob(directory + "/ana_ux_t=*.bin")[0]
    uy_file_ana: str = glob.glob(directory + "/ana_uy_t=*.bin")[0]
    ux_ana: arrayXXd = load_eigen_matrix(ux_file_ana)
    uy_ana: arrayXXd = load_eigen_matrix(uy_file_ana)
    t = extract_time_from_filename(ux_file_ana)
    plt.figure()
    mag_ana = np.sqrt(ux_ana**2, uy_ana**2)
    max_mag = np.max([max_mag, np.max(mag_ana)])
    plt.pcolor(mag_ana)
    plt.title(f"analytical, time = {t}/{max_time}")
    plt.clim(0, max_mag)
    plt.colorbar()
    plt.draw()

    plt.figure()
    mag_ana = np.sqrt((ux_ana - ux) ** 2, (uy_ana - uy) ** 2)
    max_mag = np.max([max_mag, np.max(mag_ana)])
    plt.pcolor(mag_ana)
    plt.title(f"error = {t}/{max_time}")
    plt.colorbar()
    plt.show()
    input("press Enter to exit")
    return 0


if __name__ == "__main__":
    err: int = main()
    exit(err)
