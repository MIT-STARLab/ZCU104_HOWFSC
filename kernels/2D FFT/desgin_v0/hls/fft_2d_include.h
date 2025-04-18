#ifndef _H_FFT_2D_KRNL_INCLUDE_H_
#define _H_FFT_2D_KRNL_INCLUDE_H_


// includes 
#include <hls_stream.h>
#include <hls_fft.h>
#include <complex>


using namespace std;


#define MAT_ROWS 1024
#define MAT_COLS 1024

#define MAT_SIZE (MAT_ROWS * MAT_COLS)



/********   FFT STATIC CONFIGUTATION PARAMETERS   *********/
const unsigned FFT_INPUT_WIDTH = 32;                             // sizeof(float)
const unsigned FFT_OUTPUT_WIDTH = FFT_INPUT_WIDTH;               // sizeof(float)

const unsigned FFT_ROWS_NFFT_MAX = 10;                           // log2(MAT_ROWS)
const unsigned FFT_COLS_NFFT_MAX = 10;                           // log2(MAT_COLS)

const bool     HAS_NFFT = false;                                 // has runtime configurable length? No

const unsigned FFT_PHASE_FACTOR_WIDTH =  24;                     // internal phase factor precision (24-25 for floats)
const unsigned ORDERING_OPT = hls::ip_fft::natural_order;        // Output data order (natural or bit-reversed?) bit-reversed is default

const unsigned FFT_ROWS_STAGES_BLOCK_RAM = 1;                    // (max_nfft < 10) ? 0 : (max_nfft - 9) Defines the number of block RAM stages used in the implementation.
const unsigned FFT_COLS_STAGES_BLOCK_RAM = 1;                    // (max_nfft < 10) ? 0 : (max_nfft - 9) Defines the number of block RAM stages used in the implementation.



/***************       Data Types          ************/
typedef float data_t;

typedef complex<data_t>  cmpx_data_t;



struct apra_2dfft_config_row : hls::ip_fft::params_t {
    static const unsigned input_width        = FFT_INPUT_WIDTH;
    static const unsigned output_width       = FFT_OUTPUT_WIDTH;
    // static const unsigned config_width = 8;
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
    // static const unsigned config_width = 8;
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
// void data_mover(cmpx_data_t *input_mat,
//                 cmpx_data_t *output_mat,
//                 cmpx_data_t *temp_mat,
//                 hls::stream<cmpx_data_t> &matrix_to_row_strm,
//                 hls::stream<cmpx_data_t> &matrix_from_row_strm, 
//                 hls::stream<cmpx_data_t> &matrix_to_col_strm,
//                 hls::stream<cmpx_data_t> &matrix_from_col_strm);


void stream_to_fft_row(cmpx_data_t *input_mat, hls::stream<cmpx_data_t> &matrix_to_row_strm);
void stream_from_fft_row_to_fft_col(cmpx_data_t *temp_mat, hls::stream<cmpx_data_t> &matrix_from_row_strm,  hls::stream<cmpx_data_t> &matrix_to_col_strm);
void stream_from_fft_col(cmpx_data_t *output_mat, hls::stream<cmpx_data_t> &matrix_from_col_strm);


// Row Wise FFT
void fft_row_top(bool direction, hls::stream<cmpx_data_t> &matrix_to_row_strm, hls::stream<cmpx_data_t> &matrix_from_row_strm);

void fft_row_unit(bool direction, cmpx_data_t row_in[MAT_COLS], cmpx_data_t row_out[MAT_COLS]);
void fft_row_init(bool direction, configRow_t* config_row);
void read_row_in(hls::stream<cmpx_data_t> &matrix_to_row_strm, cmpx_data_t row_in[MAT_COLS]);
void write_row_out(hls::stream<cmpx_data_t> &matrix_from_row_strm, cmpx_data_t row_out[MAT_COLS]);


// Col Wise FFT 
void fft_col_top(bool direction, hls::stream<cmpx_data_t> &matrix_to_col_strm, hls::stream<cmpx_data_t> &matrix_from_col_strm);

void fft_col_unit(bool direction, cmpx_data_t col_in[MAT_ROWS], cmpx_data_t col_out[MAT_ROWS]);
void fft_col_init(bool direction, configCol_t* config_col);
void read_col_in(hls::stream<cmpx_data_t> &matrix_to_col_strm,cmpx_data_t col_in[MAT_ROWS]);
void write_col_out(hls::stream<cmpx_data_t> &matrix_from_col_strm, cmpx_data_t col_out[MAT_ROWS]);


// Top Level 2D FFT Kernel
extern "C" void fft_2d(bool direction, cmpx_data_t *input_mat, cmpx_data_t *temp_mat, cmpx_data_t *output_mat);


#endif