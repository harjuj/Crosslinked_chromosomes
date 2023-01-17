cmake_minimum_required(VERSION 3.6)
project(crosslinked_chromosome)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -lstdc++fs")

include_directories(/home/janni/Code/boost_1_79_0)
include_directories(/home/janni/Code/eigen-3.4.0)

set(SOURCE_FILES main.cpp initialization.h moves.h global.h RandomGenerator.h crosslinks.h)
add_executable(Forward_crosslinks ${SOURCE_FILES})
