# py_action.py

from logging import NullHandler
import sys
import time
from typing import List, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from enum import Enum
from multiprocessing import shared_memory
from multiprocessing.shared_memory import SharedMemory
import ctypes as cs

from numpy._typing import NDArray

shared_memory_name: str = "py_action_memory"


class PythonAction(Enum):
    SHUTDOWN = -1
    PRINT = 0
    FIGURE = 1
    CLOSE_FIGURE = 2
    PCOLOR = 3


def pipe_readline(semph: NDArray) -> str:
    print(f"python readline semph: {semph}")
    read: str = sys.stdin.readline().strip()
    py_semaphore_give(semph)
    return read


def py_semaphore_read(semph: NDArray) -> Tuple[int, int]:
    if semph[0] < 0:
        print(f"ERROR <plot.py>: Python action semaphore has negative value")
        return (0, -9)
    return (semph[0], 0)


def py_semaphore_give(semph: NDArray) -> int:
    semph_value, err = py_semaphore_read(semph)
    if err < 0:
        print(f"ERROR <plot.py>: Failed to read python action semaphore")
        return err
    semph[0] = semph_value + 1
    return 0


def py_semaphore_take(semph: NDArray, max_wait_ms: int) -> int:
    semph_value, err = py_semaphore_read(semph)
    if err < 0:
        print(f"ERROR <plot.py>: Failed to read python action semaphore")
        return err
    start_time = time.time()
    while semph_value < 1:
        time.sleep(0.001)
        semph_value, err = py_semaphore_read(semph)
        if err < 0:
            print(f"ERROR <plot.py>: Failed to read python action semaphore")
            return err
        if max_wait_ms != -1:
            elapsed_time = (time.time() - start_time) * 1000
            if elapsed_time >= max_wait_ms:
                print(
                    f"ERROR <plot.py>: Python action semaphore take wait time exceeded"
                )
                return -10
    semph[0] = semph_value - 1
    return 0


def read_shared_data(shm: SharedMemory, Nx: int, Ny: int, Nz: int) -> NDArray:
    data: NDArray
    if Nz == 1:
        if Ny == 1:
            data = np.ndarray((Nx,), offset=8, dtype=np.float64, buffer=shm.buf)
        else:
            # INFO : the order of Nx and Ny is not a mistake (Zeile zerst, spalte später)
            data = np.ndarray(
                (Ny, Nx), offset=8, dtype=np.float64, buffer=shm.buf, order="F"
            )
    else:
        data = np.ndarray(
            (Nx, Ny, Nz), offset=8, dtype=np.float64, buffer=shm.buf, order="F"
        )
    return data


def print_action(shm: SharedMemory, semph: NDArray) -> int:
    err: int = py_semaphore_take(semph, -1)
    if err < 0:
        print(f"ERROR <plot.py>: Failed to take semaphore in print action")
        return err
    Nx: int
    Ny: int
    Nz: int
    try:
        read: str = pipe_readline(semph)
        Nx, Ny, Nz = map(int, read.split())
    except Exception as e:
        print(f"ERROR <plot.py>: Print action failed: {e}", file=sys.stderr)
        return -4
    data: NDArray = read_shared_data(shm, Nx, Ny, Nz)
    print(data)
    err = py_semaphore_give(semph)
    if err < 0:
        print(f"ERROR <plot.py>: Failed to give semaphore in print action")
        return err
    return 0


def close_figure_action(fig_idx: int) -> int:
    try:
        plt.close(fig_idx)
    except Exception as e:
        print(
            f"ERROR <plot.py>: Failed to close figure with idx {fig_idx}: {e}",
            file=sys.stderr,
        )
        return -5
    return 0


def figure_action(semph: NDArray) -> int:
    fig_idx: int
    width: float
    height: float
    figsize: Tuple[float, float]
    try:
        read: str = pipe_readline(semph)
        fig_idx = int(read)
        read = pipe_readline(semph)
        width, height = map(float, read.split())
        figsize = (width, height)
    except Exception as e:
        print(f"ERROR <plot.py>: figure action failed: {e}", file=sys.stderr)
        return -6
    plt.figure(fig_idx, figsize=figsize)
    return 0


def pcolor_action(shm: SharedMemory, semph: NDArray) -> int:
    err: int = py_semaphore_take(semph, -1)
    if err < 0:
        print(f"ERROR <plot.py>: Failed to take semaphore in pcolor action")
        return err
    fig_idx: int
    Nx: int
    Ny: int
    try:
        read: str = pipe_readline(semph)
        fig_idx = int(read)
        read = pipe_readline(semph)
        Nx, Ny = map(int, read.split())
    except Exception as e:
        print(f"ERROR <plot.py>: pcolor action failed: {e}", file=sys.stderr)
        return -7
    plt.figure(fig_idx)
    data: NDArray = read_shared_data(shm, Nx, Ny * 3, 1)
    Z: NDArray = data[:Ny, :]
    X: NDArray = data[Ny : Ny * 2, :]
    Y: NDArray = data[Ny * 2 :, :]
    plt.pcolor(X, Y, Z, shading="auto", cmap="viridis")
    plt.show()

    err = py_semaphore_give(semph)
    if err < 0:
        print(f"ERROR <plot.py>: Failed to give semaphore in pcolor action")
        return err
    return 0


def main() -> int:
    try:
        start_str: str = ""
        while start_str != "start":
            start_str = str(sys.stdin.readline().strip())
    except Exception as e:
        print(f"Error <plot.py>: Failed to start: {e}")
        return -1
    shm: SharedMemory
    shared_memory_flag: bool = False
    while shared_memory_flag is False:
        try:
            shm = shared_memory.SharedMemory(name=shared_memory_name)
        except FileNotFoundError:
            # Shared memory not available yet, wait a bit and retry
            time.sleep(0.1)
            continue
        shared_memory_flag = True
    semph: NDArray = np.ndarray((1,), dtype=np.int64, buffer=shm.buf)
    shutdown: bool = False
    while not shutdown:
        try:
            print(f"before read")
            read = pipe_readline(semph)
            print(f"Read: {read}")
            time.sleep(1)
            try:
                action: int = int(read)
            except:
                continue
        except Exception as e:
            print(f"Error <plot.py>: Failed to get action: {e}")
            return -2
        match action:
            case PythonAction.SHUTDOWN.value:
                shutdown = True
            case PythonAction.FIGURE.value:
                err = figure_action(semph)
                if err < 0:
                    return err
            case PythonAction.PRINT.value:
                err = print_action(shm, semph)
                if err < 0:
                    return err
            case PythonAction.PCOLOR.value:
                err = pcolor_action(shm, semph)
                if err < 0:
                    return err
            case _:
                print(f"ERROR <plot.py>: Python Action {action} not implemented")
                return -3
    return 0


if __name__ == "__main__":
    status: int = main()
    exit(status)
