# plot.py

import sys
import time
from typing import List, Tuple
import numpy as np
import matplotlib.pyplot as plt
from enum import Enum
from multiprocessing import shared_memory

from numpy._typing import NDArray

shared_memory_name: str = "py_action_memory"


class PythonAction(Enum):
    SHUTDOWN = -1
    PRINT = 0
    PLOT = 1
    QUIVER = 2
    IMAGESC = 3


def print_action(shm) -> int:
    try:
        Nx, Ny, Nz = map(int, sys.stdin.readline().strip().split())
    except Exception as e:
        print(f"ERROR <plot.py>: Print action failed: {e}", file=sys.stderr)
        return -4
    data: NDArray = np.ndarray((Nx, Ny, Nz), dtype=np.float64, buffer=shm.buf)
    print(data)
    return 0


def main() -> int:
    try:
        start_str: str = ""
        while start_str != "start":
            start_str = str(sys.stdin.readline().strip())
    except Exception as e:
        print(f"Error <plot.py>: Failed to start: {e}")
        return -1
    shm = None
    while shm is None:
        try:
            shm = shared_memory.SharedMemory(name=shared_memory_name)
        except FileNotFoundError:
            # Shared memory not available yet, wait a bit and retry
            time.sleep(0.1)
        shutdown: bool = False
        while not shutdown:
            try:
                read = sys.stdin.readline().strip()
                print(f"action: {read}")
                if read == "":
                    print("got empty string")
                    continue
                action: int = int(read)
            except Exception as e:
                print(f"Error <plot.py>: Failed to get action: {e}")
                return -2
            match action:
                case PythonAction.SHUTDOWN.value:
                    shutdown = True
                case PythonAction.PRINT.value:
                    err = print_action(shm)
                    if err < 0:
                        return err
                case _:
                    print(f"ERROR <plot.py>: Python Action {action} not implemented")
                    return -3
    return 0


if __name__ == "__main__":
    status: int = main()
    exit(status)
