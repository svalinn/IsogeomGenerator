cmake_minimum_required(VERSION 3.1)

# project name
project(nested_vols)


include_directories("/home/kkiesling/Pokeball/install/moab-v5.0/include")

# use c++11 standards
set (CMAKE_CXX_STANDARD 11)

add_executable(closed_nested closed_nested.cpp)
