#include "py_action.h"
#include "helpers.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <errhandlingapi.h>
#include <handleapi.h>
#include <iostream>
#include <memory>
#include <string>
#include <thread>

int py_semph_take(std::shared_ptr<int> semph_ptr, const int &max_wait_ms);
int py_semph_give(std::shared_ptr<int> semph_ptr);

#ifdef __linux__

int createSharedMemory(unique_ptr_shared_memory &sharedmemory_struct,
                       const std::size_t &shared_memory_size,
                       const char *shared_memory_name) {
  int shm_fd = shm_open(shared_memory_name, O_CREAT | O_RDWR, 0666);
  if (shm_fd == -1) {
    print_err("py_action.cpp", "Could not create shared memory (linux)");
    return 1;
  }

  // Set the size of the shared memory
  if (ftruncate(shm_fd, shared_memory_size) == -1) {
    print_err("py_action.cpp", "Error setting size of shared memory (linux)");
    close(shm_fd);
    return 1;
  }

  // Map the shared memory into the process's address space
  void *map = mmap(nullptr, shared_memory_size, PROT_READ | PROT_WRITE,
                   MAP_SHARED, shm_fd, 0);
  if (map == MAP_FAILED) {
    print_err("py_action", "Error mapping shared memory");
    close(shm_fd);
    return 1;
  }
  sharedmemory_struct->map = map;
  sharedmemory_struct->mapFile = shm_fd;
  sharedmemory_struct->size = shared_memory_size;
  sharedmemory_struct->name = shared_memory_name;
  return 0;
}

#elif _WIN32

int createSharedMemory(unique_ptr_shared_memory &sharedmemory_struct,
                       const std::size_t &shared_memory_size,
                       const char *shared_memory_name) {
  HANDLE hMapFile =
      CreateFileMappingA(INVALID_HANDLE_VALUE, nullptr, PAGE_READWRITE, 0,
                         shared_memory_size, shared_memory_name);

  if (hMapFile == NULL) {
    print_err("py_action.cpp",
              "Could not create file mapping object (windows)");
    return 1;
  }

  void *map =
      MapViewOfFile(hMapFile, FILE_MAP_ALL_ACCESS, 0, 0, shared_memory_size);

  if (map == nullptr) {
    print_err("py_action.cpp", "Could not map view of file (windows)");
    CloseHandle(hMapFile);
    return 1;
  }
  sharedmemory_struct->map = map;
  sharedmemory_struct->mapFile = hMapFile;
  sharedmemory_struct->size = shared_memory_size;
  return 0;
}

#else
print_err(
    "py_action.cpp",
    "Shared memory mapping currently only supported on _WIN32 and __linux__");
#endif

// ########### Main Program #########

const char *shared_memory_name = "py_action_memory";

int read_semph(int &semph_value, const int *semph_ptr) {
  const int read_value = *semph_ptr;
  semph_value = read_value;
  if (semph_value < 0) {
    print_err("py_action.cpp", "Python action semaphore has a negative value");
    return 1;
  }
  return 0;
}

void py_semph_init(std::shared_ptr<int> semph_ptr, const int &init_val) {
  *semph_ptr = init_val;
}

int py_semph_give(std::shared_ptr<int> semph_ptr) {
  int semph_value;
  int err = read_semph(semph_value, semph_ptr.get());
  if (err != 0) {
    print_err("py_action.cpp", "Failed to read python action semaphore");
    return 1;
  }
  *semph_ptr = semph_value + 1;
  return 0;
}

int py_semph_take(std::shared_ptr<int> semph_ptr, const int &max_wait_ms) {
  int semph_value;
  int err = read_semph(semph_value, semph_ptr.get());
  if (err != 0) {
    print_err("py_action.cpp", "Failed to read python action semaphore");
    return 1;
  }
  auto start_time = std::chrono::steady_clock::now();
  while (semph_value < 1) {
    if (max_wait_ms != -1) {
      auto current_time = std::chrono::steady_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(
          current_time - start_time);
      if (elapsed_time.count() >= max_wait_ms) {
        print_err("py_action.cpp",
                  "Python action semaphore take wait time exceeded");
        return 1;
      }
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    err = read_semph(semph_value, semph_ptr.get());
    if (err != 0) {
      print_err("py_action.cpp", "Failed to read python action semaphore");
      return 1;
    }
  }
  *semph_ptr = semph_value - 1;
  return 0;
}

int py_init(unique_ptr_shared_memory &sharedmemory_struct,
            std::shared_ptr<double[]> shared_dArray,
            std::shared_ptr<int> semph_ptr, FILE *&pipe,
            const std::size_t &mmap_max_vector_size) {

  std::string python_cmd =
      std::string("python ") + std::string(PYTHON_DIR) + "/py_action.py";

  pipe = _popen(python_cmd.c_str(), "w");
  if (!pipe) {
    shared_dArray.reset();
    std::cerr
        << "ERROR <py_action.cpp>: Failed to start Python script (command:\""
        << python_cmd << "\")" << std::endl;
    return 1;
  }
  const std::size_t shared_memory_size =
      (mmap_max_vector_size + 1) * sizeof(double);

  createSharedMemory(sharedmemory_struct, shared_memory_size,
                     shared_memory_name);

  char *base_ptr = static_cast<char *>(sharedmemory_struct->map);
  shared_dArray.reset(reinterpret_cast<double *>(base_ptr) + sizeof(int));
  semph_ptr.reset(static_cast<int *>(sharedmemory_struct->map));
  // py_semph_init(semph_ptr, 1);

  // std::size_t wait_seconds = std::max(mmap_max_vector_size / (5 * 1024 *
  // 1024),
  //                                     static_cast<std::size_t>(2));
  // fflush(stdout);
  // fflush(pipe);
  // std::this_thread::sleep_for(std::chrono::seconds(wait_seconds));
  // fprintf(pipe, "start\n");
  // fflush(pipe);
  // py_semph_give(semph_ptr);
  std::cout << "\n end of py_init \n" << std::endl;

  return 0;
}

int py_clean(unique_ptr_shared_memory &sharedmemory_struct, FILE *&pipe,
             std::shared_ptr<double[]> shared_dArray,
             std::shared_ptr<int> semph_ptr) {

  py_semph_take(semph_ptr, -1);
  fprintf(pipe, "%i\n", SHUTDOWN);
  fflush(pipe);
  py_semph_take(semph_ptr, -1);
  py_semph_give(semph_ptr);

  int status = _pclose(pipe);
  if (status < 0) {
    std::cerr << "ERROR <py_action.cpp>: Failed to close Python script pipe. "
              << "Errorcode: " << status << std::endl;
    return 1;
  }
  shared_dArray.reset();
  semph_ptr.reset();
  sharedmemory_struct.reset();
  return 0;
}

int py_print(FILE *&pipe, std::shared_ptr<double[]> shared_dArray,
             std::shared_ptr<int> semph_ptr, const Gridsize &gridsize,
             const Eigen::ArrayXXd &Z) {

  std::size_t size = gridsize[0] * gridsize[1] * gridsize[2];
  Eigen::VectorXd Z_flat = Z.reshaped();

  py_semph_take(semph_ptr, -1);
  std::memcpy(shared_dArray.get(), Z.data(), size * sizeof(double));
  py_semph_give(semph_ptr);

  py_semph_take(semph_ptr, -1);
  fprintf(pipe, "%i\n", PRINT);
  fflush(pipe);
  py_semph_take(semph_ptr, -1);
  py_semph_give(semph_ptr);

  py_semph_take(semph_ptr, -1);
  fprintf(pipe, "%zu %zu %zu\n", gridsize[0], gridsize[1], gridsize[2]);
  fflush(pipe);
  py_semph_take(semph_ptr, -1);
  py_semph_give(semph_ptr);

  return 0;
}

int py_figure(FILE *&pipe, const std::size_t &fig_idx,
              const std::array<double, 2> figsize) {
  fprintf(pipe, "%i", FIGURE);
  fflush(pipe);
  fprintf(pipe, "%zu", fig_idx);
  fflush(pipe);
  fprintf(pipe, "%f %f", figsize[0], figsize[1]);
  fflush(pipe);
  return 0;
}

int py_pcolor(FILE *&pipe, std::shared_ptr<double[]> shared_dArray,
              std::shared_ptr<int> semph_ptr, std::size_t &fig_idx,
              const Gridsize &gridsize, const Grid &grid,
              const Eigen::ArrayXXd &Z) {

  if (gridsize[2] != 1) {
    print_err("py_action.cpp",
              "py_print Function only works for 1D or 2D arrays");
    return 1;
  }
  std::size_t size = gridsize[0] * gridsize[1];
  py_semph_take(semph_ptr, -1);
  std::memcpy(shared_dArray.get(), Z.data(), size * sizeof(double));
  std::memcpy(shared_dArray.get() + size, grid.X.data(), size * sizeof(double));
  std::memcpy(shared_dArray.get() + size * 2, grid.Y.data(),
              size * sizeof(double));
  py_semph_give(semph_ptr);

  py_semph_take(semph_ptr, -1);
  fprintf(pipe, "%i\n", PCOLOR);
  fflush(pipe);
  py_semph_take(semph_ptr, -1);
  py_semph_give(semph_ptr);

  py_semph_take(semph_ptr, -1);
  fprintf(pipe, "%zu\n", fig_idx);
  fflush(pipe);
  py_semph_take(semph_ptr, -1);
  py_semph_give(semph_ptr);

  py_semph_take(semph_ptr, -1);
  fprintf(pipe, "%zu %zu\n", gridsize[0], gridsize[1]);
  fflush(pipe);
  py_semph_take(semph_ptr, -1);
  py_semph_give(semph_ptr);

  return 0;
}
