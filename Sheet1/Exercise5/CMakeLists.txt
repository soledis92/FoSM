cmake_minimum_required(VERSION 3.6)
project(exercise5_cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(exercise5_cpp ${SOURCE_FILES})