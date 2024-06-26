
/*
 * MIT STAR Lab
 * Modified by Subhi in Jun 24
 */


#include "hls_fft.h"
#include <complex>


using namespace std;


// Configurable parameters (Static)
const char FFT_INPUT_WIDTH = 32;                             // Use 32 for presecion / 34 also available  (don't fully understand)
const char FFT_OUTPUT_WIDTH = FFT_INPUT_WIDTH;

const char FFT_PHASE_FACTOR_WIDTH = 24;                      // (don't fully understand)
const char FFT_STAGES_BLOCK_RAM = 2;                         // From UG: (max_nfft < 10) ? 0 : (max_nfft - 9)
const char FFT_CONFIG_WIDTH = 16;                            // width of these params (don't fully understand)

const char FFT_NFFT_MAX = 11;                                // log2(2048)
const int FFT_LENGTH = 1 << FFT_NFFT_MAX;                    // 2048



struct my_fft_config:  hls::ip_fft::params_t {
    static const unsigned input_width   = FFT_INPUT_WIDTH;
    static const unsigned output_width  = FFT_OUTPUT_WIDTH;
    static const unsigned max_nfft      = FFT_NFFT_MAX;
    static const unsigned phase_factor_width = FFT_PHASE_FACTOR_WIDTH;
    static const unsigned stages_block_ram = FFT_STAGES_BLOCK_RAM;
    static const unsigned hax_nfft      = false;                             //since we know the specifc size ahead
    static const unsigned ordering_opt  = hls::ip_fft::natural_order;
    static const unsigned config_width  = FFT_CONFIG_WIDTH;
};


// FFT Runtime Configuration and Status
typedef hls::ip_fft::config_t<my_fft_config> config_t;                        
typedef hls::ip_fft::status_t<my_fft_config> status_t;


// Data data_type
typedef complex<float> cmpx_data_t;



/**
* Transforms the input data from array to stram to be passed to the fft module.
* @param direction   1 indicates FFT and 0 indicated InverseFFT
* @param config_s    the stream to hold the run time configuration of the fft kernel
* @param in          input data arrray
* @param out_s       output stream to be passed as an input stream to the FFT
*/
void input_data_mover(bool direction, hls::stream<config_t> &config_s, cmpx_data_t *in, hls::stream<cmpx_data_t> &out_s);



/**
* Transforms the stream output from fft module to an array when ready.
* @param status_in_s   the stream to hold the computation status
* @param ovflo         ovflo indicator
* @param in            
* @param out_s       
*/
void output_data_mover(hls::stream<status_t>  &status_in_s, bool* ovflo, hls::stream<cmpx_data_t> &in_s, cmpx_data_t *out);



 /**
  * Takes a stream of data and retruns a stream of the FFT output
  * @param direction   1 indicates FFT and 0 indicated InverseFFT
  * @param ovflo       If any overflow occurred during the FFT
  * @param in          
  * @param out_s 
  */
extern "C" void fft_top(bool direction, cmpx_data_t in[FFT_LENGTH], cmpx_data_t OUT[FFT_LENGTH], bool* ovflo);
