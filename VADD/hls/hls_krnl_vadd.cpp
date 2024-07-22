/*
# Copyright © 2023 Advanced Micro Devices, Inc. All rights reserved.
# SPDX-License-Identifier: MIT
*/

#include "vadd_include.h"


// TRIPCOUNT identifier
const int c_size = DATA_SIZE;


static void load_input(const data_t* in, hls::stream<data_t>& inStream, int size) {

    mem_rd:
        for (int i = 0; i < size; i++) {
            #pragma HLS LOOP_TRIPCOUNT min = c_size max = c_size
            #pragma HLS PIPELINE
            inStream << in[i];
        }

}



static void compute_add(hls::stream<data_t>& in1_stream,
                        hls::stream<data_t>& in2_stream,
                        hls::stream<data_t>& out_stream,
                        int size) {

    execute:
        for (int i = 0; i < size; i++) {
            #pragma HLS LOOP_TRIPCOUNT min = c_size max = c_size
            #pragma HLS PIPELINE
            out_stream << (in1_stream.read() + in2_stream.read());
        }

}



static void store_result(data_t* out, hls::stream<data_t>& out_stream, int size) {

    mem_wr:
        for (int i = 0; i < size; i++) {
            #pragma HLS LOOP_TRIPCOUNT min = c_size max = c_size
            #pragma HLS PIPELINE
            out[i] = out_stream.read();
        }

}


/*
    Vector Addition Kernel

    Arguments:
        in1  (input)  --> Input vector 1
        in2  (input)  --> Input vector 2
        out  (output) --> Output vector
        size (input)  --> Number of elements in vector
*/
void krnl_vadd(const data_t* in1, const data_t* in2, data_t* out, int size) {

    static hls::stream<data_t> in1_stream("input_stream_1");
    static hls::stream<data_t> in2_stream("input_stream_2");
    static hls::stream<data_t> out_stream("output_stream");

    #pragma HLS INTERFACE mode=m_axi bundle=gmem0 depth=c_size port=in1
    #pragma HLS INTERFACE mode=m_axi bundle=gmem1 depth=c_size port=in2
    #pragma HLS INTERFACE mode=m_axi bundle=gmem0 depth=c_size port=out
    #pragma HLS dataflow

    // dataflow pragma instruct compiler to run following three APIs in parallel
    load_input(in2, in2_stream, size);
    load_input(in1, in1_stream, size);
    compute_add(in1_stream, in2_stream, out_stream, size);
    store_result(out, out_stream, size);
}
