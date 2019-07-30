cmake_minimum_required(VERSION 3.14)
project(bioalg_cpp)

set(CMAKE_CXX_STANDARD 17)

include_directories("/home/fanta/.local/boost_1_70_0/include")
add_executable(bioalg_cpp chapter05/align.cpp)
# add_executable(bioalg_cpp chapter05/scratch.cpp)