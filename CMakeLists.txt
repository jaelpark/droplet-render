cmake_minimum_required(VERSION 2.8)
project(droplet-render)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_CXX_FLAGS "-std=c++11 -mcmodel=medium -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(Python_ADDITIONAL_VERSIONS 3.5 3.6)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

if(USE_AVX2)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")
	add_definitions(-DUSE_AVX2)
	add_definitions(-DUSE_SSE4)
	MESSAGE(STATUS "Enabling AVX2 support")
elseif(USE_SSE4)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4")
	add_definitions(-DUSE_SSE4)
	MESSAGE(STATUS "Enabling SSE4.2 support")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
	MESSAGE(STATUS "Enabling SSE3 support (default)")
endif()

if(USE_EMBREE)
	set(OPT_LIBS ${OPT_LIBS} embree)
	add_definitions(-DUSE_EMBREE)
	MESSAGE(STATUS "Enabling embree occlusion testing")
endif()

add_definitions(-DDROPLET_CPU)

file(GLOB SOURCES src/*.cpp)

add_library(droplet SHARED ${SOURCES})
TARGET_LINK_LIBRARIES(droplet tbb tbbmalloc openvdb ${OPT_LIBS})
