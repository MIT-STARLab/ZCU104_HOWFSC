#ifndef _MATRIXVECTOR_H_
#define _MATRIXVECTOR_H_
#include <hls_vector.h>

//streaming is about 4x slower in simulation
#define STREAMING 0
// Matrix is M x N
// make them all multiples of 8, and unroll loop by 8
#define MATRIX_ROWS 34*4 //(34*34)    //(ouput states DM 34*34 -> 128*128) //1158, must be a multiple of RAKEWIDTH
#define MATRIX_COLS 16*2 //(64*64*2)     //8192 (input states sensor 64*64*2 -> 200*200*2) (sensor outputs phase and amplitude)
#define MATRIX_RAKE 8 //(34) //chunk size -- how many to parallelize, must be a factor of MATRIX_ROWS
// 1156x8192 DM is 34x34, and pupil image is 64x64

typedef double data_t; //ap_int<20>

extern "C" {
//avoid memory contention
#if STREAMING
#define MEMORY_DWIDTH 512
#define MEMORY_DWIDTH_BYTES (MEMORY_DWIDTH/8)
#define NUM_WORDS 8 //((MEMORY_DWIDTH) / (8 * sizeof(data_t)))  // 64 bytes worth of doubles at a time, 8 bits/byte
#define ALLOC_SIZE(x) (NUM_WORDS*(1+x/NUM_WORDS))
// typedef struct {data_t d[NUM_WORDS];} data_chunk; // move 512 bits at a time
void matrixvector(data_t *C,const data_t *A,const data_t *B);
#else
#define ALLOC_SIZE(x) x
void matrixvector(data_t C[MATRIX_ROWS],const data_t A[MATRIX_ROWS*MATRIX_COLS],const data_t B[MATRIX_COLS]);
#endif                         
}
#endif

