
#include "vadd_new.h"
#include "hls_print.h"


void vadd_new(
            const unsigned int *in1, // Read-Only Vector 1
            const unsigned int *in2, // Read-Only Vector 2
            unsigned int *out       // Output Result
            )
{
        #pragma HLS INTERFACE m_axi port=in1 bundle=aximm1 depth=28
        #pragma HLS INTERFACE m_axi port=in2 bundle=aximm2 depth=28
        #pragma HLS INTERFACE m_axi port=out bundle=aximm1 depth=28

        for(int i = 0; i < 7; ++i)
        {
            out[i] = in1[i] + in2[i];
            hls::print("elememnt 1 is %d\n", in1[i]);
            hls::print("elememnt 2 is %d\n", in2[i]);
            hls::print("output is %d\n", out[i]);

        }
    }
