cmake_minimum_required(VERSION 3.23)
project(P1)

set(CMAKE_CXX_STANDARD 17)

#Prerequisites
find_package(Eigen3 REQUIRED)

add_executable(P1 main.cpp molecule.cpp molecule.h mass.h hf.cpp hf.h general.cpp general.h)
target_link_libraries(P1 Eigen3::Eigen)