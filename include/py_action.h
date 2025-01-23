#ifndef PY_ACTION_H
#define PY_ACTION_H

#include <cstddef>
#include <cstdio>
#include <memory>

typedef enum { PRINT, PLOT, QUIVER, IMAGESC } PythonAction;

#ifdef __linux__
typedef struct {
  void *map;
  int mapFile;
  std::size_t size;
  const char *name;
} SharedMemory_struct;

#elif _WIN32
#include <windows.h>

typedef struct {
  void *map;
  HANDLE mapFile;
} SharedMemory_struct;
#endif

int py_init(SharedMemory_struct &sharedmemory_struct,
            std::unique_ptr<double[]> &shared_dArray, FILE *&pipe,
            const std::size_t &mmap_max_vector_size);

int py_print(FILE *&pipe, std::unique_ptr<double[]> &shared_dArray,
             const std::size_t &Nx, const std::size_t &Ny,
             const std::size_t &Nz);

int py_clean(SharedMemory_struct &sharedmemory_struct, FILE *&pipe,
             std::unique_ptr<double[]> &shared_dArray);

#endif // PY_ACTION_H
