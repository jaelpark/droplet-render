#ifndef MAIN_H
#define MAIN_H

__attribute__ ((aligned(16)))
#include <xmmintrin.h>
#include <smmintrin.h> //SSE4
#include <immintrin.h> //AVX2
#include <math.h> //powf
#include "smmath.inl"
#include "SMMathPort.inl"

#include <vector>
#include <atomic>

#include <tbb/tbb.h>

typedef unsigned int uint;

//#define DebugPrintf printf
void DebugPrintf(const char *, ...);

#endif
