cmake_minimum_required(VERSION 3.6)
project(CPP_Simulation)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp json11.cpp json11.hpp Helpers.cpp Helpers.h Simulation.cpp Simulation.h)
add_executable(CPP_Simulation ${SOURCE_FILES})