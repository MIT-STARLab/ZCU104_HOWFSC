#include "qr_dcmp.h"
#include <iostream>
#include "hls_math.h"

#include <stdio.h>


//For more performance: stream in vector types that fill the entire bus width


void load_j_rows_from_stream(
    stream_t& Q_in,
    stream_t& R_in,
    vector_t Qj_ram[],
    vector_t Rj_ram[],
    const int n_vector){
    LOAD_J_ROWS_FROM_STREAM:
    for (int k=0; k<n_vector; k++){
        #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
        Qj_ram[k] = Q_in.read();
        Rj_ram[k] = R_in.read();
    }
}

void operate_on_row(
    vector_t Q_temp[],
    vector_t R_temp[],
    vector_t Qi_ram[],
    vector_t Ri_ram[],
    vector_t Qj_ram[],
    vector_t Rj_ram[],
    const int n_vector,
    const data_t sin_val,
    const data_t cos_val
){

    OPERATE_ON_ROW:
    for(int k=0; k<n_vector; k++){ //SIMD here with vectors
        #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
        #pragma HLS PIPELINE
        //#pragma HLS UNROLL factor=4
        //This can't be unrolled because each iteration depends on the last!
        R_temp[k] = Ri_ram[k];
        Ri_ram[k] = R_temp[k] * cos_val + Rj_ram[k] * (-sin_val);
        Rj_ram[k] = R_temp[k] * sin_val + Rj_ram[k] * cos_val;

        Q_temp[k] = Qi_ram[k];
        Qi_ram[k] = Q_temp[k] * cos_val + Qj_ram[k] * (-sin_val);
        Qj_ram[k] = Q_temp[k] * sin_val + Qj_ram[k] * cos_val;
    }
}

void store_j_rows_to_stream(
    stream_t& Q_out,
    stream_t& R_out,
    vector_t Qj_ram[],
    vector_t Rj_ram[],
    const int n_vector
){
    STORE_J_ROWS_TO_STREAM:
    for (int k=0; k<n_vector; k++){
        #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
        Q_out << Qj_ram[k];
        R_out << Rj_ram[k];
    }
}

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
    const int i,
    const int n){

    const int n_vector = n/VECTOR_SIZE;


    vector_t Q_temp[MAX_N/VECTOR_SIZE];
    vector_t R_temp[MAX_N/VECTOR_SIZE];

    vector_t Qi_ram[MAX_N/VECTOR_SIZE];
    vector_t Ri_ram[MAX_N/VECTOR_SIZE];


    // Partition the second dimension for parallel access
    // #pragma HLS ARRAY_PARTITION variable=Qj_ram type=block factor=128
    // #pragma HLS ARRAY_PARTITION variable=Rj_ram type=block factor=128

    // #pragma HLS ARRAY_PARTITION variable=Qi_ram type=block factor=128
    // #pragma HLS ARRAY_PARTITION variable=Ri_ram type=block factor=128

    // #pragma HLS ARRAY_PARTITION variable=Q_temp type=block factor=128
    // #pragma HLS ARRAY_PARTITION variable=R_temp type=block factor=128



    //#pragma HLS bind_storage variable=Qi_ram type=ram_2p


    // #pragma HLS bind_storage variable=R_temp type=ram_2p impl=uram
    // #pragma HLS bind_storage variable=Ri_ram type=ram_2p impl=uram
    // #pragma HLS bind_storage variable=Rj_ram type=ram_2p impl=uram
 //   #pragma HLS bind_storage variable=Q_temp impl=bram
 //   #pragma HLS bind_storage variable=R_temp impl=bram
    
    // DUMMY_TEST:
    // for (int k=0; k<n_vector; k++){
    //     #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
    //     Q_out << Q_in.read();
    //     R_out << R_in.read();

    //     for (int j = i + 1; j < n; j++){ //j starts at 1 for the first i=0

    //         #pragma HLS LOOP_TRIPCOUNT min=1 avg=N/2 max=N-1
    //         //#pragma HLS pipeline off
    //         //#pragma HLS PIPELINE II=150

    //         LOAD_J_ROWS_FROM_STREAM:
    //         for (int k=0; k<n_vector; k++){
    //             #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
    //             Q_out << Q_in.read();
    //             R_out << R_in.read();
    //         }
    //     }

    // }

    LOAD_I_ROW_FROM_STREAM:
    for (int k=0; k<n_vector; k++){
        #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
        Qi_ram[k] = Q_in.read();
        Ri_ram[k] = R_in.read();
    }


    const int loop_bound = n - (i+1);

    data_t a=Ri_ram[i/VECTOR_SIZE][i%VECTOR_SIZE];


    //#pragma HLS HLS DATAFLOW

    //Inter-loop dependence here is only the i vectors (stored in local ram)
    //Can be pipelined, but cannot be parallelized due to this inter-loop dependency
    //This CAN be pipelined dbecause j>i always, and j always increases
    //You can be writing one Qj while you're processing the next Qj
    MAIN_INNER_LOOP:
    //for (int j = i + 1; j < n; j++){ //j starts at 1 for the first i=0
    for (int jj = 0; jj < loop_bound; jj++){
        #pragma HLS LOOP_TRIPCOUNT min=1 avg=N/2 max=N-1
        //#pragma HLS PIPELINE II=100
        #pragma HLS UNROLL factor=2
        //#pragma HLS DATAFLOW //no good because Qj Rj can only be written in one process function
        int j = jj+i+1;

        vector_t Qj_ram[MAX_N/VECTOR_SIZE];
        vector_t Rj_ram[MAX_N/VECTOR_SIZE];
        //#pragma HLS DATAFLOW

        load_j_rows_from_stream(Q_in, R_in, Qj_ram, Rj_ram,n_vector);

        // Compute Givens rotation
        //These are rows from which we must extract a particular element
        data_t b=Rj_ram[i/VECTOR_SIZE][i%VECTOR_SIZE];


        data_t hypot_val = hls::hypotf(a,b);
        data_t cos_val = a/hypot_val;
        data_t sin_val = -b/hypot_val;

        operate_on_row(Q_temp,R_temp,Qi_ram,Ri_ram,Qj_ram,Rj_ram,n_vector,sin_val,cos_val);
        store_j_rows_to_stream(Q_out,R_out,Qj_ram,Rj_ram,n_vector);
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


static void combined_stream_in(const vector_t* input_matrix, const vector_t* input_row, stream_t &in_stream, const int row, int i, int n){
    int n_vector = n/VECTOR_SIZE;

    //This didn't help anything
    // vector_t local_buffer[MAX_N/VECTOR_SIZE];
    // #pragma HLS ARRAY_PARTITION variable=local_buffer complete dim=1  // Enables parallel access


    // // Preload data to local buffer (burst read from DRAM)
    // LOAD_ROW:
    // for (int k = 0; k < n_vector; k++) {
    //     #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
    //     #pragma HLS PIPELINE II=1
    //     local_buffer[k] = input_matrix[row*n_vector + k];  // One-time global memory access
    // }
    //const vector_t* my_input_row = input_matrix+row*n_vector;

    STREAM_ROW_IN:
    for (int k = 0; k < n_vector; k++){
        #pragma HLS PIPELINE II=1
        #pragma HLS LOOP_TRIPCOUNT max=N_VECTOR
        //in_stream << input_matrix[row*n_vector+k];
        in_stream << input_row[k];
        //in_stream << local_buffer[k];
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

static void combined_stream_out(stream_t &out_stream, vector_t* output_matrix, vector_t* output_row, int row, int i, int n){
    int n_vector = n/VECTOR_SIZE;

    //vector_t* my_output_row = output_matrix + row*n_vector;

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
        output_row[k] = out_stream.read();
    }
}


void do_calc(vector_t* Q,
            vector_t* R,
            stream_t& in_stream_Q,
            stream_t& in_stream_R,
            stream_t& out_stream_Q,
            stream_t& out_stream_R,
            const int n_vector,
            const int i,
            const int n){
    #pragma HLS dataflow 
    vector_t* Q_row = Q+i*n_vector;
    vector_t* R_row = R+i*n_vector;
    combined_stream_in(Q, Q_row, in_stream_Q, i,i,n);
    combined_stream_in(R, R_row, in_stream_R, i,i,n);

    inner_loop_calc(in_stream_Q, in_stream_R, out_stream_Q, out_stream_R,i,n);

    combined_stream_out(out_stream_Q, Q, Q_row, i,i,n);
    combined_stream_out(out_stream_R, R, R_row, i,i,n);
}


// QR DCMP Kernel
void krnl_qr_dcmp(data_t *Q0, data_t *R0, int n) {

    vector_t* Q = (vector_t*) Q0;
    vector_t* R = (vector_t*) R0;


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

        const int n_vector = n/VECTOR_SIZE;

        do_calc(Q, R, in_stream_Q, in_stream_R, out_stream_Q, out_stream_R,n_vector,i,n);   

    }
}