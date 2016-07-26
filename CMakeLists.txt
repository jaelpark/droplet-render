cmake_minimum_required(VERSION 2.8)
project(droplet-render)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3 -mcmodel=medium -msse4 -Wno-unused-parameter")

#legacy
add_definitions(-DBLCLOUD_CPU)
add_definitions(-DUSE_SSE4)

add_definitions(-DDROPLET_CPU)
add_definitions(-DDROPLET_SSE4)

file(GLOB SOURCES src/*.cpp)

add_library(droplet SHARED ${SOURCES})

#LINK_DIRECTORIES
TARGET_LINK_LIBRARIES(droplet tbb openvdb)

