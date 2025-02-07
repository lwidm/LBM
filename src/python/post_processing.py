#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from helpers import *
from ex_3 import *


##
# \brief Entry point of the `post_precessing.py` file
def main() -> int:
    err: int = ex_3(OUTPUT_DIR)
    return err


if __name__ == "__main__":
    err: int = main()
    exit(err)
