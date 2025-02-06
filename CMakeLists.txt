cmake_minimum_required(VERSION 3.15)
project(LBM VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ----- Use ccache if it exists -----
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
  message(STATUS "Using ccache: ${CCACHE_PROGRAM}")
else()
  message(WARNING "ccache not found, compilation will not be cached")
endif()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_BINARY_DIR}")
add_executable(
  ${PROJECT_NAME}
  "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/lbm_core.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/analytical.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/helpers.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/export_data.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/course_ex/ex_3.cpp")

target_precompile_headers(${PROJECT_NAME} PRIVATE
                          "${CMAKE_SOURCE_DIR}/include/pch.h")

# Python Files set(PYTHON_DIR "${CMAKE_SOURCE_DIR}/src/python")
# file(TO_CMAKE_PATH "${PYTHON_DIR}" PYTHON_DIR)
# target_compile_definitions(${PROJECT_NAME} PRIVATE
# PYTHON_DIR=\"${PYTHON_DIR}\")

# ----- Compilation configuration -----

add_library(LBM_compiler_flags INTERFACE)
target_compile_features(LBM_compiler_flags INTERFACE cxx_std_17)
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,APPLECLANG,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")

# Warnings (applied to all configurations)
target_compile_options(
  LBM_compiler_flags
  INTERFACE
    "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused>>"
    "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>")

# Release Configuration
target_compile_options(
  LBM_compiler_flags
  INTERFACE
    # "$<$<AND:${gcc_like_cxx},$<CONFIG:Release>>:$<BUILD_INTERFACE:-O3;-march=native;-ffast-math;-funroll-loops;-mtune=native;-fomit-frame-pointer>>"
    "$<$<AND:${gcc_like_cxx},$<CONFIG:Release>>:$<BUILD_INTERFACE:-O3;-march=native;-funroll-loops;-mtune=native;-fomit-frame-pointer>>"
    "$<$<AND:${msvc_cxx},$<CONFIG:Release>>:$<BUILD_INTERFACE:-O2;-GL;-arch:AVX2;-fp:fast;-Oi;-Ot;-Oy;-GT;-GS;-Ob3;-Qpar;-MP>>"
)

# Debug Configuration
if(gcc_like_cxx
   AND UNIX
   AND CMAKE_BUILD_TYPE STREQUAL "Debug")
  find_library(GOOGLE_PROFILER_LIBRARY profiler)
  if(GOOGLE_PROFILER_LIBRARY)
    message(STATUS "Found Google Profiler library: ${GOOGLE_PROFILER_LIBRARY}")
    set(USE_GOOGLE_PROFILER 1)
  else()
    message(WARNING "Google Profiler library not found.")
    set(USE_GOOGLE_PROFILER 0)
  endif()
else()
  set(USE_GOOGLE_PROFILER 0)
endif()

# INFO : march=native makes the profile very different
target_compile_options(
  LBM_compiler_flags
  INTERFACE
    "$<$<AND:${gcc_like_cxx},$<CONFIG:Debug>>:$<BUILD_INTERFACE:-O1;-g;-pg;-fno-omit-frame-pointer>>"
    "$<$<AND:${msvc_cxx},$<CONFIG:Debug>>:$<BUILD_INTERFACE:-O1;-Zi;-DEBUG;-PROFILE;-arch:AVX2>>"
)
# Debug-specific linker options
target_link_options(
  LBM_compiler_flags
  INTERFACE
  "$<$<AND:${gcc_like_cxx},$<CONFIG:Debug>>:$<BUILD_INTERFACE:$<IF:${USE_GOOGLE_PROFILER},-g;-lprofiler,-gp>>>"
  "$<$<AND:${msvc_cxx},$<CONFIG:Debug>>:$<BUILD_INTERFACE:-DEBUG>>")

target_link_libraries(${PROJECT_NAME} PUBLIC LBM_compiler_flags)

# ----- File configurations -----

set(OUTPUT_DIR "${CMAKE_SOURCE_DIR}/output")

configure_file("${CMAKE_SOURCE_DIR}/include/ProjectConfig.h.in"
               "${CMAKE_BINARY_DIR}/ProjectConfig.h" @ONLY)

# ----- Python files -----
configure_file("${CMAKE_SOURCE_DIR}/src/python/post_processing.py.in"
               "${CMAKE_BINARY_DIR}/post_processing.py")
if(UNIX)
  file(CHMOD "${CMAKE_BINARY_DIR}/post_processing.py" FILE_PERMISSIONS
       OWNER_EXECUTE)
endif()

# ----- Generate and copy compile commands and build directory gitignore -----
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  message(
    WARNING
      "The MSVC compiler is being used. Note that compile_commands.json cannot be generated with this compiler."
  )
else()
  add_custom_command(
    TARGET ${PROJECT_NAME}
    PRE_BUILD
    COMMAND
      ${CMAKE_COMMAND} -E copy_if_different
      "${CMAKE_BINARY_DIR}/compile_commands.json"
      "${CMAKE_SOURCE_DIR}/compile_commands.json")
endif()
add_custom_command(
  TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E touch "${CMAKE_SOURCE_DIR}/build/.gitignore")
# COMMAND ${CMAKE_COMMAND} -E echo "*" > "${CMAKE_SOURCE_DIR}/build/.gitignore")
file(WRITE "${CMAKE_SOURCE_DIR}/build/.gitignore" "*\n")

# ----- Generate output direcotry and add a gitignore -----
add_custom_command(
  TARGET ${PROJECT_NAME}
  PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory "${OUTPUT_DIR}"
  COMMAND ${CMAKE_COMMAND} -E touch "${OUTPUT_DIR}/.gitignore")
# COMMAND ${CMAKE_COMMAND} -E echo "*" > "${OUTPUT_DIR}/.gitignore")
file(WRITE "${OUTPUT_DIR}/.gitignore" "*\n")
