cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(WirePlane LANGUAGES C CXX)

# Find packages
if (NOT TARGET ROOT)
    find_package(ROOT REQUIRED)
endif()
if (NOT TARGET Garfield)
    find_package(Garfield REQUIRED)
endif()
if (NOT TARGET OpenMP)
    find_package(OpenMP REQUIRED)
endif()

# Use GNU-style directory hierarchy
include(GNUInstallDirs)

# make install output dir
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/install")

# Set the C++ standard
set(CMAKE_CXX_FLAGS -std=c++${ROOT_CXX_STANDARD})

# List of targets
set(INSTALLABLE_TARGETS "")

# ---Build libraries-------------------------------------------------------------

add_library(WirePlane SHARED src/WirePlane.cpp)
target_include_directories(WirePlane PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(WirePlane Garfield::Garfield)

# ---Build executables-----------------------------------------------------------

add_executable(field src/field.cpp)
target_link_libraries(field Garfield::Garfield ROOT::Core WirePlane)