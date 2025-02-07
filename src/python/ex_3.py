import re
import os
from helpers import *


def extract_parameter_from_directory(dirname: str) -> float:
    if "fixed_u0" in dirname:
        match = re.search(r"nu=(\d+\.\d+e[+-]\d+)", dirname)
        if match:
            return float(match.group(1))
        else:
            raise ValueError(f"Could not extract viscosity from directory: {dirname}")
    elif "fixed_nu" in dirname:
        match = re.search(r"u0=(\d+\.\d+e[+-]\d+)", dirname)
        if match:
            return float(match.group(1))
        else:
            raise ValueError(f"Could not extract velocity from directory: {dirname}")
    else:
        raise ValueError(
            "fDirectory name does not match known parameter information: {dirname}"
        )


def ex3_plot(directory: str, name: str) -> int:
    X: arrayXXd = load_eigen_matrix(directory + "/X.bin")
    Y: arrayXXd = load_eigen_matrix(directory + "/Y.bin")

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
    mag_arrays: List = [None] * len(ux_files)
    for i in np.arange(len(ux_files)):
        ux = load_eigen_matrix(ux_files[i])
        uy = load_eigen_matrix(uy_files[i])
        mag_arrays[i] = np.sqrt(ux**2, uy**2)
        max_mag = np.max([max_mag, np.max(mag_arrays[i])])

    fig = plt.figure()
    pcolor_numpy_series(name, fig, X, Y, ux_files, mag_arrays, (0, max_mag), 0.1)

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
    plt.pcolor(X, Y, mag_ana)
    plt.title(f"{name}: analytical, time = {t}/{max_time}")
    plt.clim(0, max_mag)
    plt.colorbar()
    plt.draw()

    plt.figure()
    mag_ana = np.sqrt((ux_ana - ux) ** 2, (uy_ana - uy) ** 2)
    max_mag = np.max([max_mag, np.max(mag_ana)])
    plt.pcolor(X, Y, mag_ana)
    plt.title(f"{name}: error = {t}/{max_time}")
    plt.colorbar()
    plt.show()

    plt.pause(2)
    # input("press Enter to exit")

    return 0


def ex_3(output_dir: str) -> int:

    directories: List[str] = [
        d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))
    ]
    directories = [d for d in directories if "2D_Taylor_Green" in d]
    fixed_u0_directories: List[str] = [d for d in directories if "fixed_u0" in d]
    fixed_nu_directories: List[str] = [d for d in directories if "fixed_nu" in d]
    fixed_u0_directories = sorted(
        fixed_u0_directories, key=extract_parameter_from_directory
    )
    fixed_nu_directories = sorted(
        fixed_nu_directories, key=extract_parameter_from_directory
    )

    name: str
    err: int = 0
    for directory in fixed_u0_directories:
        directory = output_dir + "/" + directory
        name = "fixed u0, nu=" + str(extract_parameter_from_directory(directory))
        err += ex3_plot(directory, name)

    for directory in fixed_nu_directories:
        directory = output_dir + "/" + directory
        name = "fixed nu, u0=" + str(extract_parameter_from_directory(directory))
        err += ex3_plot(directory, name)

    return err
