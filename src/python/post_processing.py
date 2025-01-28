import numpy as np
import matplotlib.pyplot as plt

# from numpy.typing import NDArray
# from typing import Tuple
from enum import Enum

arrayXXXd = np.ndarray[tuple[int, int, int], np.dtype[np.float64]]
arrayXXd = np.ndarray[tuple[int, int], np.dtype[np.float64]]
arrayXd = np.ndarray[tuple[int], np.dtype[np.float64]]


class Log_level(Enum):
    LOG_LEVEL_OFF = 0
    LOG_LEVEL_ERR = 1
    LOG_LEVEL_WRN = 2
    LOG_LEVEL_INF = 3
    LOG_LEVEL_DBG = 4


# %% %%%%%%%%%%%%%%%%%%%% Definitions ####################

LOG_LEVEL: Log_level = Log_level.LOG_LEVEL_DBG

filename: str = "post_processing.py"


# %% %%%%%%%%%%%%%%%%%%%% Logging functions ####################


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


def LOG_ERR(filename: str, message: str) -> None:
    log(filename, message, Log_level.LOG_LEVEL_ERR)
    return


def LOG_WRN(filename: str, message: str) -> None:
    log(filename, message, Log_level.LOG_LEVEL_WRN)
    return


def LOG_INF(filename: str, message: str) -> None:
    log(filename, message, Log_level.LOG_LEVEL_INF)
    return


def LOG_DBG(filename: str, message: str) -> None:
    log(filename, message, Log_level.LOG_LEVEL_DBG)
    return


# %% %%%%%%%%%%%%%%%%%%%% Loading matrices from files functions ####################


def load_eigen_vector(filename: str) -> arrayXd:
    matrix: arrayXd
    with open(filename, "rb") as f:
        matrix = np.fromfile(f, dtype=np.float64)
    return matrix


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


def main() -> int:
    ux_file: str = "../../output/test2/ux_t=1.0000e+03.bin"
    ux = load_eigen_matrix(ux_file)
    uy_file: str = "../../output/test2/uy_t=1.0000e+03.bin"
    uy = load_eigen_matrix(uy_file)
    mag = np.sqrt(ux**2 + uy**2)
    plt.pcolor(mag)
    plt.show()
    return 0


if __name__ == "__main__":
    err: int = main()
    exit(err)
