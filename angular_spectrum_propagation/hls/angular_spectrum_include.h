/*
 * MIT STAR Lab
 * M.Subhi Abo Rdan (msubhi_a@mit.edu)
 * Last modified in Dec 29, 2024
 */

#ifndef _H_FFT_2D_KRNL_INCLUDE_H_
#define _H_FFT_2D_KRNL_INCLUDE_H_

// includes 
#include <hls_stream.h>
#include <hls_fft.h>
#include <hls_vector.h>
#include <hls_math.h>
#include <complex>


using namespace hls;


/********   Matrix Parameters   *********/
#define MAT_ROWS 64
#define MAT_COLS 64
#define MAT_SIZE (MAT_ROWS * MAT_COLS)




/********   FFT STATIC CONFIGUTATION PARAMETERS   *********/
const unsigned FFT_INPUT_WIDTH = 32;                             // sizeof(float)
const unsigned FFT_OUTPUT_WIDTH = FFT_INPUT_WIDTH;               // sizeof(float)

const unsigned FFT_ROWS_NFFT_MAX = 6;                           // log2(MAT_ROWS)
const unsigned FFT_COLS_NFFT_MAX = 6;                           // log2(MAT_COLS)

const unsigned CONFIG_WIDTH = 8;
const bool     HAS_NFFT = false;                                 // has runtime configurable length? No to save resources

const unsigned FFT_PHASE_FACTOR_WIDTH =  24;                     // internal phase factor precision (24-25 for floats)
const unsigned ORDERING_OPT = hls::ip_fft::natural_order;        // Output data order (natural or bit-reversed?) bit-reversed is default

const unsigned FFT_ROWS_STAGES_BLOCK_RAM = 1;                    // (max_nfft < 10) ? 0 : (max_nfft - 9) Defines the number of block RAM stages used in the implementation.
const unsigned FFT_COLS_STAGES_BLOCK_RAM = 1;                    // (max_nfft < 10) ? 0 : (max_nfft - 9) Defines the number of block RAM stages used in the implementation.

const unsigned COMPLEX_MULT_TYPE = hls::ip_fft::use_mults_performance;
const unsigned BUTTERFLY_TYPE    = hls::ip_fft::use_xtremedsp_slices;


struct apra_2dfft_config_row : hls::ip_fft::params_t {
    static const unsigned input_width        = FFT_INPUT_WIDTH;
    static const unsigned output_width       = FFT_OUTPUT_WIDTH;
    static const unsigned config_width       = CONFIG_WIDTH;
    static const unsigned max_nfft           = FFT_ROWS_NFFT_MAX;
    static const bool     has_nfft           = HAS_NFFT;
    static const unsigned phase_factor_width = FFT_PHASE_FACTOR_WIDTH;
    static const unsigned ordering_opt       = ORDERING_OPT;
    static const unsigned scaling_opt        = hls::ip_fft::block_floating_point;
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
    static const unsigned scaling_opt        = hls::ip_fft::block_floating_point;
    static const unsigned stages_block_ram   = FFT_COLS_STAGES_BLOCK_RAM;
    static const unsigned complex_mult_type  = COMPLEX_MULT_TYPE;
    static const unsigned butterfly_type     = BUTTERFLY_TYPE;
};

typedef hls::ip_fft::config_t<apra_2dfft_config_row> configRow_t;
typedef hls::ip_fft::status_t<apra_2dfft_config_row> statusRow_t;
   
typedef hls::ip_fft::config_t<apra_2dfft_config_col> configCol_t;
typedef hls::ip_fft::status_t<apra_2dfft_config_col> statusCol_t;



/***************       Data Types          ************/
#define VECTOR_SIZE_ROW 8       // vectorize data streaming from/to row-wise fft; elements come from the same row
#define VECTOR_SIZE_COL 4       // vectorize data streaming and computations in col-wise fft elements come from different cols

#define FFT_CORE_ROW_COUNT 2                    // number of cores working on row-wise fft in pipelined execution
#define FFT_CORE_COL_COUNT (VECTOR_SIZE_COL)    // number of cores working on col-wise fft in parallel  execution



typedef float data_t;
typedef complex<data_t>  cmpx_data_t;

typedef hls::vector<cmpx_data_t, VECTOR_SIZE_ROW> vector_row_data_t;
typedef hls::vector<cmpx_data_t, VECTOR_SIZE_COL> vector_col_data_t;

typedef hls::stream<cmpx_data_t> cmpx_stream_t;



/***********************       Modules      ********************/
/***************************************************************/


/**
 * @brief Top-Level Module: Performs Angular Spectrum Propagation on a 2D Wavefront Input.
 * 
 * This function propagates a 2D wavefront in free space over a specified distance using the angular 
 * spectrum method. It performs forward and inverse Fourier transforms to switch between the spatial 
 * and frequency domains, applying a transfer function to model propagation. 
 * 
 * The size of the square wavefront grid (n: number of rows and columns) must be a power of 2 for 
 * FFT compatibility.
 * 
 * The following parameters are provided by the host application:
 * 
 * @param direction     Direction of the first fft applied on wavefront; must be forward
 * @param scale         Overall fft scaling effect; must be 1/(n*n)
 * @param distance      The propagation distance in meters.
 * @param k_2           Squared wavenumber
 * @param kxy           kxy[i] = ((-(double)n / 2.0 + i) + 0.5) * delkx    for 0 <= i < n
 *                      Where delkx = 2.0 * M_PI / (pixel_scale_d * n);
 *                      Where pixel_scale_d is the physical size of each pixel in the grid in units 
 *                      of [pixels/meter]. It determines the spatial resolution of the input wavefront.
 * 
 * @param input_mat     Input array of size n * n in row-major order, represents the complex-valued 2D 
 *                      wavefront in the spatial domain. Where 
 * @param output_mat    On return, it contains the propagated wavefront in the spatial domain after propagation.
 * @param temp_mat_1    Buffer of size n * n. Needed for implemention
 * @param temp_mat_2    Buffer of size n * n. Needed for implemention
 * 
 */

extern "C" 
void angular_spectrum(
    bool direction,
    double scale,
    double distance,
    double k_2,
    double *kxy,
    cmpx_data_t *input_mat,
    cmpx_data_t *output_mat,
    cmpx_data_t *temp_mat_1, 
    cmpx_data_t *temp_mat_2
    );


/**
 * TODO: specs
 */
void propagate_wave(
    double scale,
    double distance,
    double k_2,
    double *kxy,
    hls::stream<vector_col_data_t> &first_output_matrix_col_major_strm,
    hls::stream<vector_col_data_t> &second_input_matrix_col_major_strm);

/////////      DATA FLOW CONTROL        //////////////
//////////////////////////////////////////////////////
/**
 * TODO: specs
 */
void shifted_stream_to_first_fft(cmpx_data_t *input_mat, hls::stream<vector_row_data_t> &shifted_matrix_to_first_fft_strm);
void stream_from_fft_row_to_fft_col(cmpx_data_t *temp_mat,  hls::stream<vector_row_data_t> &matrix_from_row_strm,  hls::stream<vector_col_data_t> &matrix_to_col_strm);
void stream_from_fft_col_to_fft_row(cmpx_data_t *temp_mat,  hls::stream<vector_col_data_t> &matrix_from_col_strm,  hls::stream<vector_row_data_t> &matrix_to_row_strm);
void shifted_resolve_stream_from_second_fft(cmpx_data_t *output_mat, hls::stream<vector_row_data_t> &second_output_matrix_row_major_strm);



/////////////////         FFT        /////////////////
//////////////////////////////////////////////////////

// Row Wise FFT
/**
 * TODO: specs
 */
void fft_row_top(bool direction, hls::stream<vector_row_data_t> &matrix_to_row_strm, hls::stream<vector_row_data_t> &matrix_from_row_strm);
void serialize_rows_in(hls::stream<vector_row_data_t> &matrix_to_row_strm, cmpx_stream_t &row_in1, cmpx_stream_t &row_in2);
void serialize_rows_out(hls::stream<vector_row_data_t> &matrix_from_row_strm, cmpx_stream_t &row_out1, cmpx_stream_t &row_out2);
void fft_row_unit(bool direction, cmpx_stream_t &row_in, cmpx_stream_t &row_out);
void fft_row_init(bool direction, hls::stream<configRow_t> &config_row);


// Col Wise FFT 
/**
 * TODO: specs
 */

void fft_col_top(bool direction, hls::stream<vector_col_data_t> &matrix_to_col_strm, hls::stream<vector_col_data_t> &matrix_from_col_strm);
void distibute_cols_in(hls::stream<vector_col_data_t> &matrix_to_col_strm, cmpx_stream_t &col_in1, cmpx_stream_t &col_in2, cmpx_stream_t &col_in3, cmpx_stream_t &col_in4);
void distibute_cols_out(hls::stream<vector_col_data_t> &matrix_from_col_strm, cmpx_stream_t &col_out1, cmpx_stream_t &col_out2, cmpx_stream_t &col_out3, cmpx_stream_t &col_out4);
void fft_col_unit(bool direction, cmpx_stream_t &col_in, cmpx_stream_t &col_out);
void fft_col_init(bool direction, hls::stream<configCol_t> &config_col);


#endif