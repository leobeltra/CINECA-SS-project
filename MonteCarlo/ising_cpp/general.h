#pragma once
// Define precision: 1 for double, 0 for float
#define USE_DOUBLE_PRECISION 1

#if USE_DOUBLE_PRECISION
typedef double real_t;
#else
typedef float real_t;
#endif

#define DEBUG false

#define RUN_CPU false

#ifndef NUM_EXE_TO_MEASURE
#define NUM_EXE_TO_MEASURE 10
#endif