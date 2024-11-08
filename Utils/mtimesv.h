#ifndef _MTIMESV_H_
#define _MTIMESV_H_

#include <iostream>
#include <fstream>
// #include "hls_stream.h"
// #include "ap_int.h"
using namespace std;

// Matrix is M x N
// 1156x8192 DM is 34x34, and pupil image is 64x64
#define M (32)    //(ouput states DM 34*34 -> 128*128) //1158, must be a multiple of RAKEWIDTH
#define N (128)     //8192 (input states sensor 64*64*2 -> 200*200*2) (sensor outputs phase and amplitude)
#define NUM_ITERATIONS 1
#define RAKEWIDTH 16 //chunk size -- how many to parallelize

typedef double data_t; //ap_int<20>
extern "C" {
// C = A * B, just use the default memory transfer stuff
void matrixTimesVector1(data_t C[M], data_t A[M*N], data_t B[N]);

//stored in column-major order
void matrixTimesVector2(data_t C[M], data_t A[M*N], data_t B[N]);

//avoid memory contention
void matrixTimesVector3(data_t C[M], data_t A[M*N], data_t B[N]);
}
#endif

