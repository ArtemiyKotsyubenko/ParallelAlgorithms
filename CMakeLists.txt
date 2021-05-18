cmake_minimum_required(VERSION 3.17)
project(ParallelAlgorithms)

set(CMAKE_CXX_STANDARD 11)

find_package(MPI REQUIRED)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

SET(MYCOMPILE_FLAGS "-pthread -fopenmp")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${MYCOMPILE_FLAGS}")

add_executable(ParallelAlgorithms heat_mpi.cpp)