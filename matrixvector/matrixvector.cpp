#include "matrixvector.h"
#include <hls_stream.h>

#define NUM_RAKE (MATRIX_ROWS/MATRIX_RAKE)

//column-major A, in blocks of size RAKEWIDTH
//C = A * B
const int csize = MATRIX_ROWS*MATRIX_COLS;
#if STREAMING
extern "C" void do_calc(hls::stream<data_t>& C, hls::stream<data_t>& A, hls::stream<data_t>& B)
#else
extern "C" void matrixvector(data_t C[MATRIX_ROWS], const data_t A[MATRIX_ROWS*MATRIX_COLS], const data_t B[MATRIX_COLS])
#endif
{
    unsigned i,j,k,l,h;
	data_t ram[MATRIX_COLS];
	data_t cache[MATRIX_RAKE];
	data_t acc;

#if STREAMING
#else //depth is for the simulation to know how big a fifo to have
#pragma HLS INTERFACE m_axi port=C bundle=gmem0 depth=csize
#pragma HLS INTERFACE m_axi port=A bundle=gmem1 depth=csize
#pragma HLS INTERFACE m_axi port=B bundle=gmem0 depth=csize
#pragma HLS ARRAY_PARTITION variable=A factor=8 type=cyclic 
#pragma HLS ARRAY_PARTITION variable=B factor=8 type=cyclic 
#pragma HLS ARRAY_PARTITION variable=C factor=8 type=cyclic 
    // #pragma HLS INTERFACE mode=ap_memory port=A
    // #pragma HLS INTERFACE mode=ap_memory port=B
    // #pragma HLS INTERFACE mode=ap_memory port=C
#endif
	LOOP_PRELOAD: for(j=0; j < MATRIX_COLS; j++)
    {
#pragma HLS UNROLL factor=8
#if STREAMING
    	ram[j] = B.read();
#else
        ram[j] = B[j];         
#endif
    }

	LOOP_RAKE: for(i=0,k=0,h=0; i < NUM_RAKE; i++)
	{
		LOOP_STARTUP: for(j=0; j < MATRIX_RAKE; j++) 
        {
#pragma HLS UNROLL
            cache[j] = 0.0;
        }
	
		LOOP_OUTER: for(j=0; j < MATRIX_COLS; j++) 
		{
//#pragma HLS pipeline II=1  
            data_t x = ram[j]; // helps the compiler unroll around l instead of perceiving conflicts with j
            LOOP_INNER: for (l=0; l < MATRIX_RAKE; l++, k++)
            {
#pragma HLS PIPELINE II=1
#pragma HLS UNROLL factor=8
#if STREAMING
				cache[l] += A.read() * x;
#else
				cache[l] += A[k] * x;
#endif
			}
		}
		LOOP_DONE: for (j=0; j<MATRIX_RAKE; j++, h++)
		{	
#pragma HLS PIPELINE
#if STREAMING            
            C << cache[j];
#else
			C[h] = cache[j];    
#endif
        }        
	}
}

#if STREAMING

//size means the number of doubles, and will ignore any fractional units left over
static void load_input(const data_t* in, hls::stream<data_t>& inStream, unsigned size) 
{
    unsigned i,j;
    // data_chunk x;
    // unsigned vSize = size / NUM_WORDS;
    // if (size % NUM_WORDS) vSize++;

    mem_rd:  for (i = 0; i < size; i++) 
    {
#pragma HLS UNROLL factor=8
            inStream << in[i];
            // for (j=0; j< NUM_WORDS; j++) inStream << x.d[j];
    }
}

static void store_result(data_t* out, hls::stream<data_t>& out_stream, unsigned size) 
{
    unsigned i,j;
    // unsigned vSize = size / NUM_WORDS;
    // if (size % NUM_WORDS) vSize++;
    // data_chunk x;
    mem_wr: for (i = 0; i < size; i++) 
    {
#pragma HLS UNROLL factor=8
        // for (j=0; j<NUM_WORDS && size; j++, size--) x.d[j] = out_stream.read();
        out[i] = out_stream.read();
    }
}

extern "C" void matrixvector(data_t *C, const data_t *A, const data_t *B)
{
#pragma HLS INTERFACE m_axi port=C bundle=gmem0  depth=csize
#pragma HLS INTERFACE m_axi port=A bundle=gmem1  depth=csize
#pragma HLS INTERFACE m_axi port=B bundle=gmem0  depth=csize
    static hls::stream<data_t > in1("input_stream_1"); 
    static hls::stream<data_t > in2("input_stream_2");
    static hls::stream<data_t > out("output_stream");

#pragma HLS dataflow
    load_input(A, in1, MATRIX_ROWS*MATRIX_COLS);
    load_input(B, in2, MATRIX_COLS);
    do_calc(out, in1, in2);
    store_result(C, out, MATRIX_ROWS);
}
#endif