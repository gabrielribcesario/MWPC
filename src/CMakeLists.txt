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

# ---Build executables-----------------------------------------------------------
add_compile_options("-std=c++17")

add_executable(gasfile gasfile.cpp)
target_link_libraries(gasfile Garfield::Garfield ROOT::RIO)

add_executable(field field.cpp)
target_link_libraries(field Garfield::Garfield ROOT::RIO)
target_include_directories(field PRIVATE "${CMAKE_CURRENT_LIST_DIR}/../include")

add_executable(mwpc mwpc.cpp)
target_compile_options(mwpc PRIVATE "-O3;-fopenmp")
target_include_directories(mwpc PRIVATE "${CMAKE_CURRENT_LIST_DIR}/../include")
target_link_libraries(mwpc Garfield::Garfield ROOT::RIO ROOT::Tree OpenMP::OpenMP_CXX)