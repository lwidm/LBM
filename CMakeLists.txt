cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(LBM VERSION 1.0)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_executable(
  ${PROJECT_NAME}
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/analytical.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/helpers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/py_action.cpp)

# Python Files
set(PYTHON_DIR "${CMAKE_SOURCE_DIR}/src/python")
file(TO_CMAKE_PATH "${PYTHON_DIR}" PYTHON_DIR)
target_compile_definitions(${PROJECT_NAME} PRIVATE PYTHON_DIR=\"${PYTHON_DIR}\")

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
      ${CMAKE_BINARY_DIR}/compile_commands.json
      ${CMAKE_SOURCE_DIR}/compile_commands.json)
  add_custom_command(
    TARGET ${PROJECT_NAME}
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_BINARY_DIR}/.gitignore)
endif()
