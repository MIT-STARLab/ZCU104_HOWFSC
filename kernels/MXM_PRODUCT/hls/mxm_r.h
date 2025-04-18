#ifndef _MATRIXPRODUCT_H_
#define _MATRIXPRODUCT_H_


#include <hls_vector.h>
#include <hls_stream.h>


//These only affect the test bench operations
#define MATRIX_M 1024
#define MATRIX_N 1024
#define MATRIX_P 1024


#define TILE_SIZE 32
#define VECTOR_SIZE 8 // 512 bits 32 bits per float
#define N_VECTOR TILE_SIZE/VECTOR_SIZE

typedef double data_t; //in case you want to test it on integers, or wierd stuff like ap_int<20>
typedef hls::vector<data_t, VECTOR_SIZE> vector_t;
typedef hls::stream<vector_t> stream_t;

//typedef hls::stream<data_t> stream_t; // for streaming 8 things

//start_j needs to be a multiple of Vector_SIZE
extern "C" void krnl_mxm_r(data_t *C0, data_t *A0, data_t *B0, 
                            int start_i, int start_k, int start_j,
                            int M, int N, int P);
               
                        
#endif