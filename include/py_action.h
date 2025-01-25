#ifndef PY_ACTION_H
#define PY_ACTION_H

#include "main.h"
#include <cstddef>
#include <cstdio>
#include <memory>

typedef enum {
  SHUTDOWN = -1,
  PRINT,
  FIGURE,
  CLOSE_FIGURE,
  PCOLOR
} PythonAction;

#ifdef __linux__
#include <fcntl.h>
#include <sys/mman.h>
typedef struct {
  void *map;
  int mapFile;
  std::size_t size;
  const char *name;
} SharedMemory_struct;

struct SharedMemory_deleter {
  void operator()(SharedMemory_struct *sharedmemory_struct) {
    munmap(sharedmemory_struct.map, sharedmemory_struct.size);
    close(sharedmemory_struct.mapFile);
    shm_unlink(sharedmemory_struct.name);
  }
}

#elif _WIN32
#include <windows.h>

typedef struct {
  void *map;
  HANDLE mapFile;
  std::size_t size;
} SharedMemory_struct;

struct SharedMemory_deleter {
  void operator()(SharedMemory_struct *sharedmemory_struct) {
    if (sharedmemory_struct->map) {
      UnmapViewOfFile(sharedmemory_struct->map);
    }
    if (sharedmemory_struct->mapFile) {
      CloseHandle(sharedmemory_struct->mapFile);
    }
    delete sharedmemory_struct;
  }
};

#endif
using unique_ptr_shared_memory =
    std::unique_ptr<SharedMemory_struct, SharedMemory_deleter>;

int py_init(unique_ptr_shared_memory &sharedmemory_struct,
            std::shared_ptr<double[]> shared_dArray,
            std::shared_ptr<int> semph_ptr, FILE *&pipe,
            const std::size_t &mmap_max_vector_size);

int py_print(FILE *&pipe, std::shared_ptr<double[]> shared_dArray,
             std::shared_ptr<int> semph_ptr, const Gridsize &gridsize,
             const Eigen::ArrayXXd &Z);

int py_figure(FILE *&pipe, const std::size_t &fig_idx,
              const std::array<double, 2> figsize);

int py_pcolor(FILE *&pipe, std::shared_ptr<double[]> shared_dArray,
              std::shared_ptr<int> semph_ptr, std::size_t &fig_idx,
              const Gridsize &gridsize, const Grid &grid,
              const Eigen::ArrayXXd &Z);

int py_clean(unique_ptr_shared_memory &sharedmemory_struct, FILE *&pipe,
             std::shared_ptr<double[]> shared_dArray,
             std::shared_ptr<int> semph_ptr);

#endif // PY_ACTION_H
