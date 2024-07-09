#include "matrixvector.h"
#include <thread>
 
//column-major A, in blocks of size RAKEWIDTH
//C = A * B
//16 and 8 about the same speed, 32 slightly slower 

#define MULTI_VECTOR_SIZE 8

// TRIPCOUNT identifier
const int rows_c = MATRIX_ROWS;
const int cols_vector_c = MATRIX_COLS/VECTORSIZE; // Doesn't affect operation, just synthesis reports. AND CORRECTNESS OF COSIMULATION
//axi port size (in pointer arguments, so probably vectors?)
const int asize = rows_c * cols_vector_c * 8;
const int bsize = cols_vector_c *8;
const int csize = rows_c / VECTORSIZE*8;


//Now here we are going to stream 
//Note size is number of data_t elements in a row
void do_calc(stream_t& C, stream_t& A,  stream_t& B,int rows, int cols)
{
    int vector_count = cols/VECTORSIZE;
    vector_t ram[MAX_COLS/VECTORSIZE];
    vector_t outram;
    vector_t x_vector;
    
LOOP_PRELOAD: for(unsigned j=0; j < vector_count; j++){ 
#pragma HLS LOOP_TRIPCOUNT min = cols_vector_c max = cols_vector_c
        ram[j] = B.read();
}

LOOP_MULTI_OUTPUT_VECTOR: for (unsigned k=0; k<rows/VECTORSIZE; k++){ //For outputing multiple vectors at a time
#pragma HLS LOOP_TRIPCOUNT min = rows_c/VECTORSIZE max = rows_c/VECTORSIZE //Runs number of rows over vector size

//#pragma HLS UNROLL factor=4 //this actually doesnt really help as long as loop outer is unrolled
LOOP_OUTER: for (unsigned j=0; j < VECTORSIZE; j++){ //Loops over a vectors worth of parallel rows, each iteration is one row
    //#pragma HLS PIPELINE
    #pragma HLS UNROLL factor=4

    vector_t manual_parallel[MULTI_VECTOR_SIZE]; //Experimentally determined
    //This should be MULTI_VECTOR_SIZE :)
    //#pragma HLS ARRAY_PARTITION variable=manual_parallel dim=1 type=cyclic factor=8 
    

    #pragma HLS LOOP_TRIPCOUNT min = VECTORSIZE max = VECTORSIZE //Runs VECTORSIZE times
    //#pragma HLS PIPELINE

    ZERO_LOOP: for (unsigned l=0;l<MULTI_VECTOR_SIZE;l++){
        #pragma HLS UNROLL
        manual_parallel[l]=0.0; 
    }
  
    LOOP_PRODUCT: for (unsigned l=0; l < vector_count; l+=8) { //Does the vector-wise dot product
    #pragma HLS LOOP_TRIPCOUNT min = cols_vector_c/MULTI_VECTOR_SIZE max = cols_vector_c/MULTI_VECTOR_SIZE
        
    //#pragma HLS UNROLL factor=4
    #pragma HLS PIPELINE II=8 
        LOOP_INNERMOST: for (unsigned m=0; m < MULTI_VECTOR_SIZE; m++){
            manual_parallel[m] += ram[l+m] * A.read();
            // if(l==0){
            //     manual_parallel[m] = ram[l+m] * A.read();
            // }else{
            //     manual_parallel[m] += ram[l+m] * A.read();
            // }
        }
  
        }

    data_t x_scalar = 0.0;
    //LOOP_COMBINER: for(unsigned l; l<8;l++){x_vector += manual_parallel[l]}
    x_vector = manual_parallel[0]+manual_parallel[1]+manual_parallel[2]+manual_parallel[3]+manual_parallel[4]+manual_parallel[5]+manual_parallel[6]+manual_parallel[7];
    x_scalar = x_vector[0]+x_vector[1]+x_vector[2]+x_vector[3]+x_vector[4]+x_vector[5]+x_vector[6]+x_vector[7];
    // LOOP_VECTOR: for (unsigned l=0; l<VECTORSIZE; l++)
    //     x_scalar += x_vector[l];
    // } //Computes the elementwise sum of the vector

    outram[j] = x_scalar;
}
    C << outram; 
}
}


static void load_input(const vector_t* in, stream_t& inStream, int size) {
mem_rd:
    for (int i = 0; i < size; i++) {
#pragma HLS LOOP_TRIPCOUNT min = cols_vector_c max = cols_vector_c
        inStream << in[i];
    }
}

//Duplicate of the above just for tripcount
static void load_input_mat(const vector_t* in, stream_t& inStream, int size) {
mem_rd:
    for (int i = 0; i < size; i++) {
#pragma HLS LOOP_TRIPCOUNT min = cols_vector_c*rows_c max = cols_vector_c*rows_c
        inStream << in[i];
    }
}


static void store_result(vector_t* out, stream_t& out_stream, int size) { //size is number of vector elements
mem_wr:
    for (int i = 0; i < size; i++) {
#pragma HLS LOOP_TRIPCOUNT min = rows_c/VECTORSIZE max = rows_c/VECTORSIZE
        out[i] = out_stream.read(); //Only going to work with one at a time
    }
}

// size is the number of data_t elements, not vector_t elements! it must be a multiple of VECTORSIZE
//A should point to the start of the matrix block (start of a row)
//B should point to the start of the vector
//The result C is not a full vector, just one element or PARALLELROWS elements
extern "C" void krnl_inner_product(data_t* C0, const data_t* A0, const data_t* B0, int rows, int cols)
{

    vector_t *C = (vector_t *)C0;
    const vector_t *A = (vector_t *)A0;
    const vector_t *B = (vector_t *)B0;

    static stream_t in1("input_stream_1"); 
    static stream_t in2("input_stream_2"); 
    static stream_t out("output_stream");
    // static data_t ram[MATRIX_COLS];
//The below may be better as an AXI stream, but we would need to be careful about sizes
#pragma HLS INTERFACE m_axi port=C bundle=gmem0 depth=csize
#pragma HLS INTERFACE m_axi port=A bundle=gmem1 depth=asize
#pragma HLS INTERFACE m_axi port=B bundle=gmem0 depth=bsize
#pragma HLS dataflow

    load_input(B, in2, cols/VECTORSIZE);

    load_input_mat(A, in1, rows * cols/VECTORSIZE);
    do_calc(out, in1, in2,rows,cols);
    store_result(C, out, rows/VECTORSIZE);
}
