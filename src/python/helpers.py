import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from enum import Enum
import re
import glob

# from numpy.typing import NDArray
from typing import List, Tuple

arrayXXXd = np.ndarray[tuple[int, int, int], np.dtype[np.float64]]
arrayXXd = np.ndarray[tuple[int, int], np.dtype[np.float64]]
arrayXd = np.ndarray[tuple[int], np.dtype[np.float64]]

OUTPUT_DIR = "@OUTPUT_DIR@"


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
# \brief Plot a time series of 2D numpy arrays using pcolor.
#
# This function plots a sequence of 2D numpy arrays as a time series using `pcolor` from `matplotlib`.
# Each array is plotted sequentially, with a pause between frames to create an animation-like effect.
#
# \param name (str): The name of the plot (will be displayed in the title).
# \param fig (Figure): The matplotlib figure to use for plotting (created with `plt.figure()`).
# \param X (arrayXXd): The X grid (2D array filled with the x-values).
# \param Y (arrayXXd): The Y grid (2D array filled with the y-values).
# \param array_filenames (List[str]): A sorted (ascending in time) list of filnames corresponding to the simulation data.
# \param arrays (List[arrayXXd]): A sorted (ascending in time) list arrays corresponding to the simulation data.
# \param clim_ (Tuple[float, float]): The colorbar limits for the pcolor plot (min, max).
# \param pause (float): The pause duratio (in seconds) between each frame to control the animation speed.
#
# \see extract_time_from_filename: Function to extract a floating point time from the data's filename
def pcolor_numpy_series(
    name: str,
    fig: Figure,
    X: arrayXXd,
    Y: arrayXXd,
    array_filenames: List[str],
    arrays: List[arrayXXd],
    clim_: Tuple[float, float],
    pause: float,
) -> int:
    t_end: float = extract_time_from_filename(array_filenames[-1])
    for i in np.arange(len(arrays)):
        plt.clf()
        if not plt.fignum_exists(fig.number):
            LOG_INF(this_filename, "plot closed before all states were plotted")
            return 0
        t: float = extract_time_from_filename(array_filenames[i])
        plt.pcolor(X, Y, arrays[i])
        plt.title(f"{name}: time = {t}/{t_end}")
        plt.clim(clim_[0], clim_[1])
        plt.colorbar()
        # plt.draw()
        fig.canvas.draw_idle()
        fig.canvas.flush_events()
        # plt.pause(pause)
    return 0
