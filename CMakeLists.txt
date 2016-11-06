cmake_minimum_required(VERSION 3.6)
project(CPP_Simulation)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

include_directories(/usr/local/opt/boost/include /opt/localsolver_6_0/include /Library/gurobi700/mac64/include)
link_directories(/usr/local/opt/boost/lib /opt/localsolver_6_0/bin /Library/gurobi700/mac64/lib)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -llocalsolver -lboost_system -lboost_filesystem -lgurobi_c++ -lgurobi70")

set(SOURCE_FILES main.cpp json11.cpp json11.hpp Helpers.cpp Helpers.h Simulation.cpp Simulation.h Evaluator.cpp Evaluator.h GurobiSolver.cpp GurobiSolver.h LSSolver.cpp LSSolver.h Stopwatch.cpp Stopwatch.h Matrix.h)
add_executable(CPP_Simulation ${SOURCE_FILES})