cmake_minimum_required(VERSION 3.10.2)
project(ivr_code)

include (CTest)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

include_directories("/Users/ayuseiya/research/ivr/ivr_code/MPICode")

SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER  mpicc)

add_executable(LW_code_clean_version main.cpp Detector.cpp FullSystem.cpp System.cpp
        AuxiliaryFunction.cppp util.h
        ./MPICode/ReadMPI.cpp
        ./MPICode/ConstructDetectorMPI.cpp ./MPICode/CommunicationFunctionMPI.cpp
        ./MPICode/DetectorEvolveMPI.cpp
        ./MPICode/SystemReadMPI.cpp ./MPICode/ComputeDetectorEnergyWindowMPI.cpp
        )

target_link_libraries(ivr_code stdc++fs mpi)