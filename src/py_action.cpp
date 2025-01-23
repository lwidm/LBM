#include "py_action.h"
#include "helpers.h"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <memory>
#include <minwindef.h>
#include <string>
#include <thread>

#ifdef __linux__
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

int createSharedMemory(SharedMemory_struct &sharedmemory_struct,
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
  sharedmemory_struct.map = map;
  sharedmemory_struct.mapFile = shm_fd;
  sharedmemory_struct.size = shared_memory_size;
  sharedmemory_struct.name = shared_memory_name;
  return 0;
}

int cleanSharedMemory(SharedMemory_struct &sharedmemory_struct) {
  munmap(sharedmemory_struct.map, sharedmemory_struct.size);
  close(sharedmemory_struct.mapFile);
  shm_unlink(sharedmemory_struct.name);
}

#elif _WIN32

int createSharedMemory(SharedMemory_struct &sharedmemory_struct,
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
  sharedmemory_struct.map = map;
  sharedmemory_struct.mapFile = hMapFile;
  return 0;
}

int cleanSharedMemory(SharedMemory_struct &sharedmemory_struct) {
  UnmapViewOfFile(sharedmemory_struct.map);
  CloseHandle(sharedmemory_struct.mapFile);
  return 0;
}
#else
print_err(
    "py_action.cpp",
    "Shared memory mapping currently only supported on _WIN32 and __linux__");
#endif

// ########### Main Program #########

#define PYTHON_DIR "C:/Users/lukas/Documents/04_projects/LBM/src/python"
std::string python_cmd = std::string("python ") + PYTHON_DIR + "/plot.py";
const char *shared_memory_name = "py_action_memory";

int py_init(SharedMemory_struct &sharedmemory_struct,
            std::unique_ptr<double[]> &shared_dArray, FILE *&pipe,
            const std::size_t &mmap_max_vector_size) {

  pipe = _popen(python_cmd.c_str(), "w");
  if (!pipe) {
    shared_dArray.reset();
    std::cerr
        << "ERROR <py_action.cpp>: Failed to start Python script (command: \""
        << python_cmd << "\")" << std::endl;
    return 1;
  }
  const std::size_t shared_memory_size = mmap_max_vector_size * sizeof(double);

  createSharedMemory(sharedmemory_struct, shared_memory_size,
                     shared_memory_name);

  shared_dArray.reset(static_cast<double *>(sharedmemory_struct.map));

  std::size_t wait_seconds = std::max(mmap_max_vector_size / (2 * 1024 * 1024),
                                      static_cast<std::size_t>(5));
  std::this_thread::sleep_for(std::chrono::seconds(wait_seconds));
  fflush(pipe);
  fflush(stdout);

  fprintf(pipe, "start\n");
  fflush(pipe);
  return 0;
}

int py_clean(SharedMemory_struct &sharedmemory_struct, FILE *&pipe,
             std::unique_ptr<double[]> &shared_dArray) {

  cleanSharedMemory(sharedmemory_struct);
  fprintf(pipe, "-1\n");
  fflush(pipe);

  int status = _pclose(pipe);
  if (status < 0) {
    std::cerr << "ERROR <py_action.cpp>: Failed to close Python script pipe. "
              << "Errorcode: " << status << " (command : \"" << python_cmd
              << "\")" << std::endl;
    return 1;
  }
  shared_dArray.reset();
  return 0;
}

int py_print(FILE *&pipe, std::unique_ptr<double[]> &shared_dArray,
             const std::size_t &Nx, const std::size_t &Ny,
             const std::size_t &Nz) {
  fprintf(pipe, "%i\n", PRINT);
  fflush(pipe);
  fprintf(pipe, "%zu %zu %zu\n", Nx, Ny, Nz);
  fflush(pipe);
  return 0;
}
