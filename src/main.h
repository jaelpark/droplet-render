#ifndef MAIN_H
#define MAIN_H

/*#define _M_X64
typedef unsigned int uint32_t;
#include <DirectXMath.h>
#include <DirectXPackedVector.h>
using namespace DirectX;*/

__attribute__ ((aligned(16)))
#include <xmmintrin.h>
#include <smmintrin.h> //SSE4
#include <math.h> //powf
#include "smmath.inl"
#include "smport.inl"

#define _CRT_RAND_S

#include <vector>
#include <atomic>

typedef unsigned int uint;

//#define BLCLOUD_DEBUG

#ifdef BLCLOUD_DEBUG
#pragma optimize("",off)
void DebugPrintf(const char *pmsg, ...);
#else
#define DebugPrintf printf
#endif

#endif
