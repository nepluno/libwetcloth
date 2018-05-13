# Cmake 2.6.0 is missing a few macros that FindEigen3.cmake uses
cmake_minimum_required (VERSION 2.6.2)

project (libWetCloth)
set (LIBWETCLOTH_VERSION "2018")
add_definitions (-DLIBWETCLOTH_VERSION="${LIBWETCLOTH_VERSION}")

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install")

# Initialize the build type (Release, Debug, etc)
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Debug CACHE STRING
    "Choose the type of build, options are: Debug Release."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

add_definitions (-DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")

# Add warnings to the compiler flags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused -std=c++14 -pthread")

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions (-DDEBUG)
endif (CMAKE_BUILD_TYPE MATCHES Debug)

if (CMAKE_BUILD_TYPE MATCHES Release)
  add_definitions (-DNDEBUG)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
endif (CMAKE_BUILD_TYPE MATCHES Release)

# Add directory with macros
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Add cmake file with useful functions
include (libWetClothCommon)

option (USE_MKL "Builds in support for mkl" OFF)
if (USE_MKL)
set (MKL_STATIC 1)
find_package (MKL REQUIRED)
if (MKL_FOUND)
  include_directories (${MKL_INCLUDE_DIRS})
  set (LIBWETCLOTH_LIBRARIES ${LIBWETCLOTH_LIBRARIES} ${MKL_LIBRARIES})
  add_definitions (-DEIGEN_USE_MKL_ALL)
else (MKL_FOUND)
  message (SEND_ERROR "Unable to locate MKL")
endif (MKL_FOUND)
endif (USE_MKL)

# Eigen3 library is required
find_package (Eigen3 REQUIRED)
if (EIGEN3_FOUND)
  include_directories (${EIGEN3_INCLUDE_DIR})
else (EIGEN3_FOUND)
  message (SEND_ERROR "Unable to locate Eigen3")
endif (EIGEN3_FOUND)

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include/)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include/eigen/)
add_subdirectory (libWetCloth)

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/assets ${CMAKE_CURRENT_BINARY_DIR}/libWetCloth/assets )

