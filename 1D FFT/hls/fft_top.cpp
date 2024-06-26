/*
 * MIT STAR Lab
 * Modified by Subhi in Jun 24
 */


#include "fft_top.h"
#include "hls_fft.h"



/**
 * Inherts docs.
 */
void input_data_mover(bool direction, hls::stream<config_t> &config_s, cmpx_data_t *in, hls::stream<cmpx_data_t> &out_s){

    //Load the runtime configuration to the stream that is going to the FFT module: runtime fongif includes direction
    config_t config;
    config.setDir(direction);
    // config.setSch(0x2AB);             // don't know what that is and it seems not required in case of floating points
    config_s.write(config);

    for (int i =0; i<FFT_LENGTH; i++){
        #pragma HLS pipeline II = 1 rewind
        out_s.write(in[i]);
    }

}


/**
 * Inherts docs.
 */
void output_data_mover(hls::stream<status_t>  &status_in_s, bool* ovflo, hls::stream<cmpx_data_t> &in_s, cmpx_data_t * out){

    for (int i=0; i<FFT_LENGTH; i++){
        #pragma HLS pipeline II = 1 rewind
        out[i] = (in_s.read());
    }

    status_t status_in = status_in_s.read();
    *ovflo = status_in.getOvflo() & 0x1;                 // Don't fully understand 

}



/**
 * Inherts docs.
 */
void fft_top(bool direction, cmpx_data_t *in, cmpx_data_t *out, bool* ovflo){

    #pragma HLS interface s_axilite  port = direction
    #pragma HLS interface s_axilite  port = ovflo
    #pragma HLS interface m_axi      port = in    depth = FFT_LENGTH
    #pragma HLS interface m_axi      port = out   depth = FFT_LENGTH

    #pragma HLS dataflow
    hls::stream<cmpx_data_t> input_stream;
    hls::stream<cmpx_data_t> output_stream;
    hls::stream<config_t> runtime_config;
    hls::stream<status_t> runtime_status;

    // convert inputs to hls::stream<> and generates fft_config stream based on input arguments
    input_data_mover(direction, runtime_config, in, input_stream);

    // FFT IP Instantiate
    hls::fft<my_fft_config>(input_stream, output_stream, runtime_status, runtime_config);

    //Move output stream to the output array
    output_data_mover(runtime_status, ovflo, output_stream, out);

}

