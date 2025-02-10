#include "qr_dcmp.h"
#include <iostream>
#include "hls_math.h"

#include <stdio.h>


//For more performance: stream in vector types that fill the entire bus width


void inner_loop_calc(
    // stream_t& Qi_in, 
    // stream_t& Qj_in, 
    // stream_t& Ri_in, 
    // stream_t& Rj_in, 

    // stream_t& Qi_out, 
    // stream_t& Qj_out, 
    // stream_t& Ri_out, 
    // stream_t& Rj_out, 
    stream_t& Q_in,
    stream_t& R_in,
    stream_t& Q_out,
    stream_t& R_out,
    int i,
    int n){

    int n_vector = n/VECTOR_SIZE;


    vector_t Q_temp[MAX_N/VECTOR_SIZE];
    vector_t R_temp[MAX_N/VECTOR_SIZE];

    vector_t Qi_ram[MAX_N/VECTOR_SIZE];
    vector_t Qj_ram[MAX_N/VECTOR_SIZE];
    vector_t Ri_ram[MAX_N/VECTOR_SIZE];
    vector_t Rj_ram[MAX_N/VECTOR_SIZE];


    // #pragma HLS bind_storage variable=R_temp type=ram_2p impl=uram
    // #pragma HLS bind_storage variable=Ri_ram type=ram_2p impl=uram
    // #pragma HLS bind_storage variable=Rj_ram type=ram_2p impl=uram
 //   #pragma HLS bind_storage variable=Q_temp impl=bram
 //   #pragma HLS bind_storage variable=R_temp impl=bram


    LOAD_I_ROW_FROM_STREAM:
    for (int k=0; k<n_vector; k++){
        #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
        Qi_ram[k] = Q_in.read();
        Ri_ram[k] = R_in.read();
    }

    data_t a;
    data_t b;
    data_t hypot_val;
    data_t cos_val;
    data_t sin_val;

    //Inter-loop dependence here is only the i vectors (stored in local ram)
    //Can be pipelined, but cannot be parallelized due to this inter-loop dependency
    MAIN_INNER_LOOP:
    for (int j = i + 1; j < n; j++){ //j starts at 1 for the first i=0

        #pragma HLS LOOP_TRIPCOUNT min=1 avg=N/2 max=N-1
        //#pragma HLS pipeline off
        //#pragma HLS PIPELINE II=150

        LOAD_J_ROWS_FROM_STREAM:
        for (int k=0; k<n_vector; k++){
            #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
            Qj_ram[k] = Q_in.read();
            Rj_ram[k] = R_in.read();
        }

        // Compute Givens rotation
        //These are rows from which we must extract a particular element
        a=Ri_ram[i/VECTOR_SIZE][i%VECTOR_SIZE];
        b=Rj_ram[i/VECTOR_SIZE][i%VECTOR_SIZE];

        hypot_val = hypotf(a,b);
        cos_val = a/hypot_val;
        sin_val = -b/hypot_val;

        #pragma HLS DATAFLOW


        OPERATE_ON_ROW:
        for(int k=0; k<n_vector; k++){ //SIMD here with vectors
            #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
            //#pragma HLS UNROLL factor=4
            //This can't be unrolled because each iteration depends on the last!
            R_temp[k] = Ri_ram[k];
            Ri_ram[k] = R_temp[k] * cos_val + Rj_ram[k] * (-sin_val);
            Rj_ram[k] = R_temp[k] * sin_val + Rj_ram[k] * cos_val;

            Q_temp[k] = Qi_ram[k];
            Qi_ram[k] = Q_temp[k] * cos_val + Qj_ram[k] * (-sin_val);
            Qj_ram[k] = Q_temp[k] * sin_val + Qj_ram[k] * cos_val;
        }

        STORE_J_ROWS_TO_STREAM:
        for (int k=0; k<n_vector; k++){
            #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
            Q_out << Qj_ram[k];
            R_out << Rj_ram[k];
        }
    }

    STORE_I_ROW_TO_STREAM:
    for (int k=0; k<n_vector; k++){
        #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
        Q_out << Qi_ram[k];
        R_out << Ri_ram[k];
    }

}
////////////////////////////     Data Moving    /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


void combined_stream_in(vector_t* input_matrix, stream_t &in_stream, int row, int i, int n){
    int n_vector = n/VECTOR_SIZE;

    STREAM_ROW_IN:
    for (int k = 0; k < n_vector; k++){
        #pragma HLS PIPELINE II=1
        #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
        in_stream << input_matrix[row*n_vector+k];
    }

    STREAM_TRIANGLE_IN:
    for (int j = i + 1; j < n; j++){
        #pragma HLS LOOP_TRIPCOUNT min=1 avg=N/2 max=N-1
        INNER_STREAM_TRIANGLE_IN:
        for (int k=0; k<n_vector; k++){
            #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
            #pragma HLS PIPELINE II=1
            in_stream << input_matrix[j*n_vector+k]; //Jth row, in row major order
        }
    }
}

void combined_stream_out(stream_t &out_stream, vector_t* output_matrix, int row, int i, int n){
    int n_vector = n/VECTOR_SIZE;

    STREAM_TRIANGLE_OUT:
    for (int j = i + 1; j < n; j++){
        #pragma HLS LOOP_TRIPCOUNT min=1 avg=N/2 max=N-1
        INNER_STREAM_TRIANGLE_OUT:
        for (int k=0; k<n_vector; k++){
            #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
            #pragma HLS PIPELINE II=1
            output_matrix[j*n_vector+k] = out_stream.read();  //Jth row, in row major order
        }
    }

    STREAM_ROW_OUT:
    for (int k = 0; k < n_vector; k++){
        #pragma HLS PIPELINE II=1
        #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
        output_matrix[row*n_vector+k] = out_stream.read();
    }
}


// QR DCMP Kernel
void krnl_qr_dcmp(data_t *Q0, data_t *R0, int n) {

    vector_t* Q = (vector_t*) Q0;
    vector_t* R = (vector_t*) R0;

    int n_vector = n/VECTOR_SIZE;

    #pragma HLS INTERFACE m_axi port=Q bundle=gmem0 depth=N*N_VECTOR
    #pragma HLS INTERFACE m_axi port=R bundle=gmem1 depth=N*N_VECTOR

    

    //Each bundle can only have one reader and one writer.
    //The i streaming is the row stream
    //The j streaming is the triangle stream    
    //We combine the row and triangle streaming operations by
    // streaming in the i row first, and stream out the i row last

    //i ranges from 0 to N-2 inclusive
    //j ranges from i+1 to N-1 inclusive

    //int i =1;
    for(int i = 0; i < n-1; i++){
        static stream_t in_stream_Q("in_stream_Q");
        static stream_t in_stream_R("in_stream_R");

        static stream_t out_stream_Q("out_stream_Q");
        static stream_t out_stream_R("out_stream_R");

       #pragma HLS LOOP_TRIPCOUNT min=N-1 max=N-1
       #pragma HLS DEPENDENCE variable=Q inter true distance=1
       #pragma HLS DEPENDENCE variable=R inter true distance=1  

        #pragma HLS dataflow 
        //stream_row_in(Q,in_stream_Qi,i);
        //stream_row_in(R,in_stream_Ri,i);

        //triangle_to_stream(Q, in_stream_Qj, i);
        //triangle_to_stream(R, in_stream_Rj, i);

        combined_stream_in(Q, in_stream_Q, i,i,n);
        combined_stream_in(R, in_stream_R, i,i,n);

        // inner_loop_calc(in_stream_Qi, in_stream_Qj, in_stream_Ri, in_stream_Rj,
        //              out_stream_Qi, out_stream_Qj, out_stream_Ri, out_stream_Rj, i);
        inner_loop_calc(in_stream_Q, in_stream_R, out_stream_Q, out_stream_R,i,n);    

        // stream_to_triangle(out_stream_Qj,Q,i);
        // stream_to_triangle(out_stream_Rj,R,i);

        // stream_row_out(out_stream_Ri, R, i);
        // stream_row_out(out_stream_Qi, Q, i);

         combined_stream_out(out_stream_Q, Q, i,i,n);
         combined_stream_out(out_stream_R, R, i,i,n);
    }
}

