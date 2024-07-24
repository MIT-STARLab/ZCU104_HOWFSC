/*
 * MIT STAR Lab
 * Modified by Subhi in Jun 24
 */

#include "fft_krnl_include.h"




/**
 * Convert inputs to hls::stream<>
 * Generates fft_config stream based on input arguments (direction)
 */
static void input_data_mover(const bool direction, hls::stream<config_t> &config_s, const cmpx_input_data_t *in, hls::stream<cmpx_input_data_t> &in_s){

    //////  Load the runtime configuration to the configuration stream that is going to the FFT module; includes direction, scalings
    config_t config;
    config.setDir(direction);
    // config.setSch(0x2AB);             // not required in case of unscaled FFT
    config_s << config;


    for (int i =0; i<DATA_SIZE; i++){
        #pragma HLS pipeline
        in_s << in[i];
    }

}


/**
 * Convert output to array
 * Gets runtime status struct
 */
static void output_data_mover(hls::stream<status_t>  &runtime_status_stream, bool ovflo, hls::stream<cmpx_output_data_t> &out_s, cmpx_output_data_t *out){
    for (int i=0; i<DATA_SIZE; i++){
        #pragma HLS pipeline
        out[i] = (out_s.read());
    }

    status_t status_in = runtime_status_stream.read();
    ovflo = status_in.getOvflo() & 0x1;    
    // ovflo = false;
}




/**
 * Takes an array of data and retruns a DFFT 
 * @param direction   1 indicates FFT and 0 indicated InverseFFT
 * @param in          input data array
 * @param out_s       output data array
 * @param ovflo       If any overflow occurred during the FFT
*/
void fft_top(const bool direction, const cmpx_input_data_t *in, cmpx_output_data_t *out, bool ovflo) {

    #pragma HLS interface s_axilite  port = direction
    #pragma HLS interface m_axi      port = in          bundle=gmem0  depth = DATA_SIZE
    #pragma HLS interface m_axi      port = out         bundle=gmem0  depth = DATA_SIZE
    #pragma HLS interface s_axilite  port = ovflo
    
    hls::stream<cmpx_input_data_t> input_stream;
    hls::stream<cmpx_output_data_t> output_stream;
    hls::stream<config_t> runtime_config_stream;
    hls::stream<status_t> runtime_status_stream;

    #pragma HLS dataflow
    // convert inputs to hls::stream<> and generates fft_config stream based on input arguments
    input_data_mover(direction, runtime_config_stream, in, input_stream);

    // FFT IP Instantiate
    hls::fft<star_lab_fft1d_config>(input_stream, output_stream, runtime_status_stream, runtime_config_stream);

    //Move output stream to the output array, and wire runtime status
    output_data_mover(runtime_status_stream, ovflo, output_stream, out);
}
