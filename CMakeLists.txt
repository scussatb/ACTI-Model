
cmake_minimum_required(VERSION 2.8.9)
project(Tests)

set(CMAKE_CXX_FLAGS "-L/usr/local/lib -I/usr/local/include -O3 -std=c++11 -Wall -Wextra -pedantic")
add_subdirectory(src)
