cmake_minimum_required(VERSION 3.6)
project(CPP_Simulation)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(GTest REQUIRED)

include_directories(${Boost_INCLUDE_DIRS} ${GTEST_INCLUDE_DIRS} /opt/localsolver_6_0/include /Library/gurobi700/mac64/include)

link_directories(/opt/localsolver_6_0/bin /Library/gurobi700/mac64/lib)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -llocalsolver -lgurobi_c++ -lgurobi70")

set(COMMON_SOURCE_FILES json11.cpp json11.hpp Helpers.cpp Helpers.h Simulation.cpp Simulation.h Evaluator.cpp Evaluator.h GurobiSolver.cpp GurobiSolver.h LSSolver.cpp LSSolver.h Stopwatch.cpp Stopwatch.h Matrix.h PSSolver.h PSSolver.cpp Runner.cpp Runner.h Experiments.cpp Experiments.h Globals.h)

set(SOURCE_FILES main.cpp ${COMMON_SOURCE_FILES})
add_executable(CPP_Simulation ${SOURCE_FILES})
target_link_libraries(CPP_Simulation ${Boost_LIBRARIES})

set(TEST_SOURCE_FILES UnitTests/testmain.cpp UnitTests/MatrixTest.h UnitTests/MatrixTest.cpp UnitTests/TestHelpers.h UnitTests/TestHelpers.cpp)
add_executable(CPP_Simulation_Test ${TEST_SOURCE_FILES})
target_link_libraries(CPP_Simulation_Test ${Boost_LIBRARIES} ${GTEST_BOTH_LIBRARIES})