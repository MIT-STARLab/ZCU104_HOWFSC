#ifndef _H_FFT_2D_KRNL_INCLUDE_H_
#define _H_FFT_2D_KRNL_INCLUDE_H_


// includes 
#include <hls_stream.h>
#include <hls_fft.h>
#include <hls_vector.h>
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

const unsigned CONFIG_WIDTH = 16;
const bool     HAS_NFFT = false;                                 // has runtime configurable length? No

const unsigned FFT_PHASE_FACTOR_WIDTH =  24;                     // internal phase factor precision (24-25 for floats)
const unsigned ORDERING_OPT = hls::ip_fft::natural_order;        // Output data order (natural or bit-reversed?) bit-reversed is default

const unsigned FFT_ROWS_STAGES_BLOCK_RAM = 1;                    // (max_nfft < 10) ? 0 : (max_nfft - 9) Defines the number of block RAM stages used in the implementation.
const unsigned FFT_COLS_STAGES_BLOCK_RAM = 1;                    // (max_nfft < 10) ? 0 : (max_nfft - 9) Defines the number of block RAM stages used in the implementation.

const unsigned COMPLEX_MULT_TYPE = hls::ip_fft::use_mults_performance;
const unsigned BUTTERFLY_TYPE = hls::ip_fft::use_xtremedsp_slices;



/***************       Data Types          ************/
#define VECTOR_SIZE_ROW 8       // when loading data to do fft on row_wise, load 8 element at a time from the row
#define VECTOR_SIZE_COL 4       // when loading data to do fft on col-wise, load 4 elements form a row, each correspond to a differnt col

#define FFT_CORE_ROW_COUNT 2
#define FFT_CORE_COL_COUNT (VECTOR_SIZE_COL)

struct apra_2dfft_config_row : hls::ip_fft::params_t {
    static const unsigned input_width        = FFT_INPUT_WIDTH;
    static const unsigned output_width       = FFT_OUTPUT_WIDTH;
    static const unsigned config_width       = CONFIG_WIDTH;
    static const unsigned max_nfft           = FFT_ROWS_NFFT_MAX;
    static const bool     has_nfft           = HAS_NFFT;
    static const unsigned phase_factor_width = FFT_PHASE_FACTOR_WIDTH;
    static const unsigned ordering_opt       = ORDERING_OPT;
    static const unsigned scaling_opt = hls::ip_fft::block_floating_point;
    static const unsigned stages_block_ram   = FFT_ROWS_STAGES_BLOCK_RAM;
    static const unsigned complex_mult_type  = COMPLEX_MULT_TYPE;
    static const unsigned butterfly_type     = BUTTERFLY_TYPE;
};

struct apra_2dfft_config_col : hls::ip_fft::params_t {
    static const unsigned input_width        = FFT_INPUT_WIDTH;
    static const unsigned output_width       = FFT_OUTPUT_WIDTH;
    static const unsigned config_width       = CONFIG_WIDTH;
    static const unsigned max_nfft           = FFT_COLS_NFFT_MAX;
    static const bool     has_nfft           = HAS_NFFT;
    static const unsigned phase_factor_width = FFT_PHASE_FACTOR_WIDTH;
    static const unsigned ordering_opt       = ORDERING_OPT;
    static const unsigned scaling_opt = hls::ip_fft::block_floating_point;
    static const unsigned stages_block_ram   = FFT_COLS_STAGES_BLOCK_RAM;
    static const unsigned complex_mult_type  = COMPLEX_MULT_TYPE;
    static const unsigned butterfly_type     = BUTTERFLY_TYPE;
};

typedef hls::ip_fft::config_t<apra_2dfft_config_row> configRow_t;
typedef hls::ip_fft::status_t<apra_2dfft_config_row> statusRow_t;
   

typedef hls::ip_fft::config_t<apra_2dfft_config_col> configCol_t;
typedef hls::ip_fft::status_t<apra_2dfft_config_col> statusCol_t;

typedef float data_t;

typedef complex<data_t>  cmpx_data_t;

typedef hls::vector<cmpx_data_t, VECTOR_SIZE_ROW> vector_row_data_t;
typedef hls::vector<cmpx_data_t, VECTOR_SIZE_COL> vector_col_data_t;

typedef hls::stream<cmpx_data_t> cmpx_stream_t;




/*****   Modules  *****/

// Data movers
void stream_to_fft_row(cmpx_data_t *input_mat, hls::stream<vector_row_data_t> &matrix_to_row_strm);
void stream_from_fft_row_to_fft_col(cmpx_data_t *temp_mat, hls::stream<vector_row_data_t> &matrix_from_row_strm,  hls::stream<vector_col_data_t> &matrix_to_col_strm);
void stream_from_fft_col(cmpx_data_t *output_mat, hls::stream<vector_col_data_t> &matrix_from_col_strm);


// Row Wise FFT
void fft_row_top(bool direction, hls::stream<vector_row_data_t> &matrix_to_row_strm, hls::stream<vector_row_data_t> &matrix_from_row_strm);

void serialize_rows_in(hls::stream<vector_row_data_t> &matrix_to_row_strm, cmpx_stream_t &row_in1, cmpx_stream_t &row_in2);
void serialize_rows_out(hls::stream<vector_row_data_t> &matrix_from_row_strm, cmpx_stream_t &row_out1, cmpx_stream_t &row_out2);

void fft_row_unit(bool direction, cmpx_stream_t &row_in, cmpx_stream_t &row_out);
void fft_row_init(bool direction, hls::stream<configRow_t> &config_row);




// Col Wise FFT 
void fft_col_top(bool direction, hls::stream<vector_col_data_t> &matrix_to_col_strm, hls::stream<vector_col_data_t> &matrix_from_col_strm);

void distibute_cols_in(hls::stream<vector_col_data_t> &matrix_to_col_strm, cmpx_stream_t &col_in1, cmpx_stream_t &col_in2, cmpx_stream_t &col_in3, cmpx_stream_t &col_in4);
void distibute_cols_out(hls::stream<vector_col_data_t> &matrix_from_col_strm, cmpx_stream_t &col_out1, cmpx_stream_t &col_out2, cmpx_stream_t &col_out3, cmpx_stream_t &col_out4);

void fft_col_unit(bool direction, cmpx_stream_t &col_in, cmpx_stream_t &col_out);
void fft_col_init(bool direction, hls::stream<configCol_t> &config_col);


// Top Level 2D FFT Kernel
extern "C" void fft_2d(bool direction, cmpx_data_t *input_mat, cmpx_data_t *temp_mat, cmpx_data_t *output_mat);


#endif