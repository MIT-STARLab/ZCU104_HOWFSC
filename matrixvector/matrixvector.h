#ifndef _MATRIXVECTOR_H_
#define _MATRIXVECTOR_H_

#include <iostream>
#include <fstream>

// Matrix is M x N
// 1156x8192 DM is 34x34, and pupil image is 64x64
#define MATRIX_ROWS (34*34)    //(ouput states DM 34*34 -> 128*128) //1158, must be a multiple of RAKEWIDTH
#define MATRIX_COLS (64*64*2)     //8192 (input states sensor 64*64*2 -> 200*200*2) (sensor outputs phase and amplitude)
#define MATRIX_RAKE (34) //chunk size -- how many to parallelize, must be a factor of MATRIX_ROWS
#define NUM_ITERATIONS 1

typedef double data_t; //ap_int<20>
extern "C" {
//avoid memory contention
void matrixvector(data_t C[MATRIX_ROWS], data_t A[MATRIX_ROWS*MATRIX_COLS], data_t B[MATRIX_COLS]);
}
#endif

