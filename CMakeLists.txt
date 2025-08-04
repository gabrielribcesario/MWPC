cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(MWPC LANGUAGES C CXX)

if (NOT TARGET ROOT)
    find_package(ROOT REQUIRED COMPONENTS RIO Tree)
endif()
if (NOT TARGET Garfield::Garfield)
    find_package(Garfield REQUIRED)
endif()
if (NOT TARGET OpenMP::OpenMP_CXX)
    find_package(OpenMP REQUIRED)
endif()

add_compile_options("-std=c++17")

include_directories("${CMAKE_SOURCE_DIR}/include" "${ROOT_INCLUDE_DIRS}")

file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build")

# ---Build executables-----------------------------------------------------------

add_executable(gasfile src/gasfile.cpp)
target_link_libraries(gasfile Garfield::Garfield ROOT::RIO)

add_executable(field src/field.cpp)
target_link_libraries(field Garfield::Garfield ROOT::RIO)

add_executable(mwpc src/mwpc.cpp)
target_compile_options(mwpc PRIVATE -O3 -fopenmp)
target_link_libraries(mwpc Garfield::Garfield ROOT::RIO ROOT::Tree OpenMP::OpenMP_CXX)