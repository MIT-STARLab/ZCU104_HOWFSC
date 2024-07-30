#ifndef _H_FFT_2D_KRNL_INCLUDE_H_
#define _H_FFT_2D_KRNL_INCLUDE_H_


// includes 
#include <hls_stream.h>
#include <hls_fft.h>
#include <complex>


using namespace std;


#define MAT_ROWS = 2048
#define MAT_COLS = 1024


/********   FFT STATIC CONFIGUTATION PARAMETERS   *********/
const unsigned FFT_INPUT_WIDTH = 32;                             // sizeof(float)
const unsigned FFT_OUTPUT_WIDTH = FFT_INPUT_WIDTH;               // sizeof(float)


const unsigned FFT_ROWS_NFFT_MAX = 11;                           // log2(MAT_ROWS)
const unsigned FFT_COLS_NFFT_MAX = 10;                           // log2(MAT_COLS)

const bool     HAS_NFFT = false;                                 // has runtime configurable length? No

const unsigned FFT_PHASE_FACTOR_WIDTH =  24;                     // internal phase factor precision (24-25 for floats)
const unsigned ORDERING_OPT = hls::ip_fft::natural_order;        // Output data order (natural or bit-reversed?) bit-reversed is default

const unsigned FFT_ROWS_STAGES_BLOCK_RAM = 2;                    // (max_nfft < 10) ? 0 : (max_nfft - 9) Defines the number of block RAM stages used in the implementation.
const unsigned FFT_COLS_STAGES_BLOCK_RAM = 1;                    // (max_nfft < 10) ? 0 : (max_nfft - 9) Defines the number of block RAM stages used in the implementation.



/***************       Data Types          ************/
typedef float data_in_t;
typedef float data_out_t;
   
typedef complex<data_in_t>  cmpxDataIn;
typedef complex<data_out_t> cmpxDataOut;



struct apra_2dfft_config_row : hls::ip_fft::params_t {
    static const unsigned input_width        = FFT_INPUT_WIDTH;
    static const unsigned output_width       = FFT_OUTPUT_WIDTH;

    static const unsigned max_nfft           = FFT_ROWS_NFFT_MAX;
    static const bool     has_nfft           = HAS_NFFT;

    static const unsigned phase_factor_width = FFT_PHASE_FACTOR_WIDTH;
    static const unsigned ordering_opt       = ORDERING_OPT;

    // static const unsigned scaling_opt = hls::ip_fft::block_floating_point;      //default
    static const unsigned stages_block_ram   = FFT_ROWS_STAGES_BLOCK_RAM;
};

typedef hls::ip_fft::config_t<apra_2dfft_config_row> configRow_t;
typedef hls::ip_fft::status_t<apra_2dfft_config_row> statusRow_t;
   

struct apra_2dfft_config_col : hls::ip_fft::params_t {
    static const unsigned input_width        = FFT_INPUT_WIDTH;
    static const unsigned output_width       = FFT_OUTPUT_WIDTH;

    static const unsigned max_nfft           = FFT_COLS_NFFT_MAX;
    static const bool     has_nfft           = HAS_NFFT;

    static const unsigned phase_factor_width = FFT_PHASE_FACTOR_WIDTH;
    static const unsigned ordering_opt       = ORDERING_OPT;

    // static const unsigned scaling_opt = hls::ip_fft::block_floating_point;      //default I guess
    static const unsigned stages_block_ram   = FFT_COLS_STAGES_BLOCK_RAM;
};


typedef hls::ip_fft::config_t<apra_2dfft_config_col> configCol_t;
typedef hls::ip_fft::status_t<apra_2dfft_config_col> statusCol_t;




/*****   Modules  *****/

// Data Moving and Control Unit
void stream_to_row_wise(cmpxDataIn *input_mat, hls::stream<cmpxDataIn> matrix_to_row);
void stream_from_row_wise_to_col_wise(cmpxDataIn *input_mat, cmpxDataOut *output_mat, hls::stream<cmpxDataOut> matrix_from_row,  hls::stream<cmpxDataIn> matrix_to_col);
void stream_from_col_wise(cmpxDataOut *output_mat, hls::stream<cmpxDataOut> matrix_from_col);


// Row Wise FFT
void row_wise_fft(bool direction, hls::stream<cmpxDataIn> matrix_to_row, hls::stream<cmpxDataOut> matrix_from_row);
void fft_row_init(bool direction, configRow_t config_row);
void read_row_in(hls::stream<cmpxDataIn> matrix_to_row, cmpxDataIn *row_in);
void write_row_out(hls::stream<cmpxDataOut> matrix_from_row, cmpxDataOut *row_out);


// Col Wise FFT 
void col_wise_fft(bool direction, hls::stream<cmpxDataIn> matrix_to_col, hls::stream<cmpxDataOut> matrix_from_col);
void fft_col_init(bool direction, configCol_t config_col);
void read_col_in(hls::stream<cmpxDataIn> matrix_to_col, cmpxDataIn *col_in);
void write_col_out(hls::stream<cmpxDataOut> matrix_from_col, cmpxDataOut *col_out);


// Top Level 2D FFT Kernel
extern "C" void fft_2d(bool direction, cmpxDataIn *input_mat, cmpxDataOut *output_mat, bool ovflo);


#endif
