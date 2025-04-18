#ifndef _H_QR_DCMP_KRNL_INCLUDE_H_
#define _H_QR_DCMP_KRNL_INCLUDE_H_


// includes 
#include <hls_stream.h>
#include <hls_vector.h>
#include <complex>


using namespace std;

#define MAX_N 2048

#define VECTOR_SIZE 16

//These two used only for synthesis reports
#define N 1024
#define N_VECTOR N/VECTOR_SIZE


#define MAT_SIZE (N * N)


typedef float data_t;




typedef hls::vector<data_t, VECTOR_SIZE> vector_t;

//typedef hls::stream<data_t> stream_t;
typedef hls::stream<vector_t> stream_t;

typedef hls::vector<data_t, N> row_t;


/*****   Modules  *****/

// Data movers
// void triangle_to_stream(data_t* input_matrix, stream_t &in_stream, int i);
// void stream_to_triangle(stream_t &out_stream, data_t* output_matrix, int i);
// void stream_row_in(data_t* input_matrix, stream_t &in_stream, int row);
// void stream_row_out(stream_t &out_stream, data_t* output_matrix, int row);

static void combined_stream_in(const vector_t* input_matrix, const vector_t* input_row, stream_t &in_stream, int row, int i, int n);
static void combined_stream_out(stream_t &out_stream, vector_t* output_matrix, vector_t* output_row, int row, int i, int n);


// Top Level QR DCMP Kernel
extern "C" void krnl_qr_dcmp(data_t *QT, data_t *R, int n);



#endif