cmake_minimum_required(VERSION 2.8)
project(summary)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
include_directories( "/usr/include/eigen3" )
set(SOURCE_FILES main.cpp)
add_executable(summary ${SOURCE_FILES})