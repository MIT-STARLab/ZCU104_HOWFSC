#ifndef _MATRIXVECTOR_H_
#define _MATRIXVECTOR_H_

// Define the matrix dimensions and vector siz
//These only affect the test bench operations
#define MATRIX_ROWS 256 // Should be divisible by PARALLELROWS
#define MATRIX_COLS 256*20 


#define VECTORSIZE 8 //how many words to transfer at once. 512 bits address bus?
//incrementing by 1 makes it easier for DATAFLOW to do block transfer on its own
#include <hls_vector.h>
#include <hls_stream.h>

#define MAX_COLS 1024*20 //in data_t


#define PARALLELROWS 8 // how many dot products to compute in kernel, must be a factor of MATRIX_ROWS and product of VECTORSIZE


typedef double data_t; //in case you want to test it on integers, or wierd stuff like ap_int<20>
typedef hls::vector<data_t,VECTORSIZE> vector_t;
typedef hls::stream<vector_t> stream_t; // for streaming 8 things

//avoid memory contention
#define MEMORY_DWIDTH 512
#define MEMORY_DWIDTH_BYTES (MEMORY_DWIDTH/8)


#define ALLOC_SIZE(x) (NUM_WORDS*(1+x/NUM_WORDS))
// typedef struct {data_t d[NUM_WORDS];} data_chunk; // move 512 bits at a time
extern "C" void krnl_inner_product(data_t *C,const data_t *A,const data_t *B, int rows, int cols);
                        
#endif