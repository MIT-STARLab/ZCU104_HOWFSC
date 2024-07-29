/*
 * Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
 * Copyright 2022-2024 Advanced Micro Devices, Inc. All Rights Reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


/*
 * MIT STAR Lab-
 * Modified by Subhi in Jul, 2024-
 */

#ifndef _H_FFT_KRNL_INCLUDE_H_
#define _H_FFT_KRNL_INCLUDE_H_


///////  Includes  ///////
#include "hls_fft.h"
#include <hls_stream.h>
#include <complex>
// #include <ap_fixed.h>



#define DATA_SIZE 2048 


/* FFT STATIC CONFIGUTATION PARAMETERS */
const unsigned FFT_INPUT_WIDTH = 32;                             // sizeof(float)
const unsigned FFT_OUTPUT_WIDTH = FFT_INPUT_WIDTH;               // sizeof(float), in case of unscaled ap_fixed: +FFT_NFFT_MAX + 1

// const unsigned CONFIG_WIDTH = 16;                                // default, might want to change incase of different arch

const unsigned FFT_NFFT_MAX = 11;                                // log2(DATA_SIZE)
const bool     HAS_NFFT = false;                                 // has runtime configurable length? No

const unsigned FFT_PHASE_FACTOR_WIDTH =  24;                     // internal phase factor precision (24-25 for floats)
const unsigned ORDERING_OPT = hls::ip_fft::natural_order;        // Output data order (natural or bit-reversed?) bit-reversed is default

// const unsigned SCALING_OPT = hls::ip_fft::unscaled;           //  for ap_fixed only; becareful with data types and port widths
const unsigned FFT_STAGES_BLOCK_RAM = 2;                         // (max_nfft < 10) ? 0 : (max_nfft - 9) Defines the number of block RAM stages used in the implementation.


struct star_lab_fft1d_config:  hls::ip_fft::params_t {
    static const unsigned input_width        = FFT_INPUT_WIDTH;
    static const unsigned output_width       = FFT_OUTPUT_WIDTH;
    // static const unsigned config_width = CONFIG_WIDTH;

    static const unsigned max_nfft           = FFT_NFFT_MAX;
    static const bool     has_nfft           = HAS_NFFT;

    static const unsigned phase_factor_width = FFT_PHASE_FACTOR_WIDTH;
    static const unsigned ordering_opt       = ORDERING_OPT;

    // static const unsigned scaling_opt        = SCALING_OPT;
    static const unsigned stages_block_ram   = FFT_STAGES_BLOCK_RAM;
};



/*
    ap_fixed, unscaled example:
    if input_width = 32
    then ouput_width = input_width + FFT_NFFT_MAX + 1

    typedef ap_fixed<FFT_INPUT_WIDTH, 1> in_data_t;     //ap_fixed<32, 1>

    // the width of output dara type must be output_width rounded to the nearest multiple of 8
    typedef ap_fixed<48, 17> out_data_t;                //ap_fixed<32 + (11 + 1) + 4 = 48, 1 + (11 + 1) + 4>
*/


// FFT Runtime Configuration and Status
typedef hls::ip_fft::config_t<star_lab_fft1d_config> config_t;                        
typedef hls::ip_fft::status_t<star_lab_fft1d_config> status_t;


typedef float in_data_t; // use different for input and output because in ap_fixed they are different
typedef float out_data_t;


typedef std::complex<out_data_t> cmpx_output_data_t;
typedef std::complex<in_data_t> cmpx_input_data_t;




/**
 * Takes an array of data and retruns a DFFT 
 * @param direction   1 indicates FFT and 0 indicated InverseFFT
 * @param in          input data array
 * @param out_s       output data array
 * @param ovflo       If any overflow occurred during the FFT
*/
extern "C" void fft_top(const bool direction, const cmpx_input_data_t *in, cmpx_output_data_t *out, bool ovflo);



#endif