cmake_minimum_required(VERSION 2.8)
project(droplet-render)

set(CMAKE_CXX_FLAGS "-std=c++11 -mcmodel=medium -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(USE_SSE4)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4")
	add_definitions(-DUSE_SSE4)
	MESSAGE(STATUS "Enabling SSE4.2 support")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
	MESSAGE(STATUS "Enabling SSE3 support (default)")
endif()

add_definitions(-DBLCLOUD_CPU)
add_definitions(-DDROPLET_CPU)

file(GLOB SOURCES src/*.cpp)

add_library(droplet SHARED ${SOURCES})
TARGET_LINK_LIBRARIES(droplet tbb openvdb)
