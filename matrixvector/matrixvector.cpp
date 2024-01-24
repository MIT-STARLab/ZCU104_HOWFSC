#include "matrixvector.h"
// #include <hls_vector.h>
// #include <hls_stream.h>

#define MEMORY_DWIDTH 512
#define NUM_RAKE (MATRIX_ROWS/MATRIX_RAKE)
// #define NUM_WORDS ((MEMORY_DWIDTH) / (8 * sizeof(data_t)))  // 64 bytes worth of doubles at a time, 8 bits/byte

//column-major A, in blocks of size RAKEWIDTH
//C = A * B
extern "C" void matrixvector(data_t C[MATRIX_ROWS], data_t A[MATRIX_ROWS*MATRIX_COLS], data_t B[MATRIX_COLS])
// void do_calc(hls::stream<data_t>& C, hls::stream<data_t>& A, hls::stream<data_t>& B)
{
    unsigned i,j,k,l,h;
	static data_t ram[MATRIX_COLS];
	static data_t cache[MATRIX_RAKE];
	data_t acc;

#pragma HLS DATAFLOW
// #pragma HLS INTERFACE mode=ap_memory port=A
// #pragma HLS INTERFACE mode=ap_memory port=B
// #pragma HLS INTERFACE mode=ap_memory port=C
#pragma HLS INTERFACE m_axi port=A bundle=gmem0
#pragma HLS INTERFACE m_axi port=B bundle=gmem1
#pragma HLS INTERFACE m_axi port=C bundle=gmem0
#pragma HLS ARRAY_PARTITION variable=A factor=32 type=cyclic 
#pragma HLS ARRAY_PARTITION variable=B factor=32 type=cyclic 
#pragma HLS ARRAY_PARTITION variable=C factor=34 type=cyclic 
	LOOP_PRELOAD: for(j=0; j < MATRIX_COLS; j++)
    {
#pragma HLS PIPELINE
         ram[j] = B[j];         
	// ram[j] = B.read();
    }

	LOOP_RAKE: for(i=0,k=0,h=0; i < NUM_RAKE; i++)
	{
		LOOP_STARTUP: for(j=0; j < MATRIX_RAKE; j++) 
        {
#pragma HLS PIPELINE
            cache[j] = 0.0;
        }
	
		LOOP_OUTER: for(j=0; j < MATRIX_COLS; j++) 
		{
#pragma HLS PIPELINE
				LOOP_INNER: for (l=0; l < MATRIX_RAKE; l++, k++)
				{
					// cache[l] += A.read() * ram[j];
#pragma HLS PIPELINE
#pragma HLS UNROLL factor=34
					cache[l] += A[k] * ram[j];
				}
		}
		LOOP_DONE: for (j=0; j<MATRIX_RAKE; j++, h++)
		{	
            // C << cache[j];
			C[h] = cache[j];    
        }        
	}
}
/*
//size means the number of doubles, and will ignore any fractional units left over
static void load_input(data_t* in, hls::stream<data_t>& inStream, int size) 
{
    hls::vector<data_t, NUM_WORDS> x ;
    int vSize = size / NUM_WORDS;

    mem_rd:  for (int i = 0; i < vSize; i++) {
    #pragma HLS LOOP_TRIPCOUNT min = c_size max = c_size
            x = in[i];
            for (int j=0; j< NUM_WORDS; j++)
                inStream << x[j]; //in[i];
        }
    //leftovers, not optimized    
    for (int i=vSize*NUM_WORDS; i < size; i++) inStream << x[i];
}

static void store_result(hls::vector<double, NUM_DOUBLES> * out, hls::stream<data_t>& out_stream, int size) 
{
 int vSize = size / NUM_WORDS;
 hls::vector<data_t, NUM_WORDS> x;
mem_wr:
    for (int i = 0; i < vSize; i++) {
#pragma HLS LOOP_TRIPCOUNT min = c_size max = c_size
    for (int j=0; j<NUM_WORDS; j++) x[j] = out_stream.read();
    out[i] = x;
    }
}

extern "C" void matrixvector(hls::vector<double, NUM_DOUBLES> *C[MATRIX_ROWS], 
                            hls::vector<double, NUM_DOUBLES> *A[MATRIX_ROWS*MATRIX_COLS], 
                            hls::vector<double, NUM_DOUBLES> *B[MATRIX_COLS])
{
#pragma HLS INTERFACE m_axi port = C bundle = gmem1
#pragma HLS INTERFACE m_axi port = A bundle = gmem0
#pragma HLS INTERFACE m_axi port = B bundle = gmem0
    static hls::stream<double > in1("input_stream_1"); 
    static hls::stream<double > in2("input_stream_2");
    static hls::stream<double > out("output_stream");

#pragma HLS dataflow
    load_input(A, in1, MATRIX_ROWS*MATRIX_COLS);
    load_input(B, in2, MATRIX_COLS);
    do_calc(out, in1, in2);
    store_result(C, out, MATRIX_ROWS);
}
*/