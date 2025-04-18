#include "mxm_r.h"
#include <iostream>



#include "mxm_r.h"
#include <iostream>



static void do_calc(stream_t& C_out_stream, stream_t& A_stream, stream_t& B_stream,  stream_t& C_stream) {
    // Vectorized on-chip buffers
    vector_t ram_A[TILE_SIZE][N_VECTOR];
    vector_t ram_B[TILE_SIZE][N_VECTOR];
    vector_t ram_C[TILE_SIZE][N_VECTOR];
    
    // Partition the second dimension for parallel access
    #pragma HLS ARRAY_PARTITION variable=ram_A dim=2 factor=4 type=block
    #pragma HLS ARRAY_PARTITION variable=ram_B dim=1 factor=4 type=block

    // #pragma HLS ALLOCATION type=operation instances=fadd limit=64
    // #pragma HLS ALLOCATION type=operation instances=fmul limit=64

    #pragma HLS bind_storage variable=ram_A type=ram_t2p
    #pragma HLS bind_storage variable=ram_B type=ram_t2p

    // #pragma HLS ARRAY_PARTITION variable=ram_A dim=2 type=complete
    // #pragma HLS ARRAY_PARTITION variable=ram_B dim=1 type=complete
    

    //-----------------------------------------
    // 1) PRELOAD: Read from streams into ram_A, ram_B, ram_C
    //-----------------------------------------
    auto preload = [&]() {
        PRELOAD:
        for (int i = 0; i < TILE_SIZE; i++) {
            #pragma HLS LOOP_TRIPCOUNT min=TILE_SIZE max=TILE_SIZE
            for (int j = 0; j < N_VECTOR; j++) {
                #pragma HLS LOOP_TRIPCOUNT min=N_VECTOR max=N_VECTOR
                #pragma HLS PIPELINE II=1
                ram_A[i][j] = A_stream.read();
                ram_B[i][j] = B_stream.read();
                ram_C[i][j] = C_stream.read();
            }
        }
    };

  
    // Now call the three local "functions" in order
    preload();

    // Outer loops: Unroll for parallelism across tile elements
    MMULT_OUTER:
    for (int i = 0; i < TILE_SIZE; i++) {
        //#pragma HLS LOOP_TRIPCOUNT min=TILE_SIZE max=TILE_SIZE
    
        MMULT_MID:
        for (int j = 0; j < N_VECTOR; j++) {
            #pragma HLS LOOP_TRIPCOUNT min=N_VECTOR max=N_VECTOR
            #pragma HLS PIPELINE II=16
                    
            vector_t temp = 0.0f;
            vector_t temp_arr[TILE_SIZE];
            
            // Inner loop over k is pipelined
            MMULT_INNER:
            for (int k = 0; k < TILE_SIZE; k++) {
                //#pragma HLS LOOP_TRIPCOUNT min=TILE_SIZE/8 max=TILE_SIZE/8
                #pragma HLS PIPELINE off
                #pragma HLS UNROLL 
                int a_vec_idx  = k / VECTOR_SIZE;
                int a_elem_idx = k % VECTOR_SIZE;
                data_t a_val   = ram_A[i][a_vec_idx][a_elem_idx];
                temp_arr[k] = a_val * ram_B[k][j];
            }
            MMULT_SUM:
            for(int k = 0; k < TILE_SIZE; k++){
                //#pragma HLS LOOP_TRIPCOUNT min=TILE_SIZE/8 max=TILE_SIZE/8
                #pragma HLS PIPELINE off
                #pragma HLS UNROLL 
                // #pragma HLS DEPENDENCE variable=temp_arr type=inter false
                temp += temp_arr[k];
            }    
            // Accumulate into C
            for (int v = 0; v < VECTOR_SIZE; v++){
                #pragma HLS PIPELINE off
                #pragma HLS UNROLL
                ram_C[i][j][v] += temp[v];
            }
            //ram_C[i][j] += temp;
            C_out_stream << ram_C[i][j];
            //C_out_stream << temp;
        }
    }

    
}


static void load_input_mat(vector_t *in, stream_t& in_stream, int start_i, int vec_start_j, int vec_cols) {
    //#pragma HLS INLINE
    MEM_RD:
    for (int i = 0; i < TILE_SIZE; i++) {
        #pragma HLS LOOP_TRIPCOUNT min=TILE_SIZE max=TILE_SIZE   
        for (int j = 0; j < N_VECTOR; j++) {
            #pragma HLS LOOP_TRIPCOUNT min=N_VECTOR max=N_VECTOR
            #pragma HLS PIPELINE II = 1
            vector_t val = in[(start_i + i) * vec_cols + (vec_start_j + j)];
            in_stream << val;
        }
    }
}


//Function to store results from stream into global memory
static void store_result(vector_t *out, stream_t& out_stream, int start_i, int vec_start_j, int vec_cols) {
    MEM_WR:
    for (int i = 0; i < TILE_SIZE; i++) {
        #pragma HLS LOOP_TRIPCOUNT min=TILE_SIZE max=TILE_SIZE      
        for (int j = 0; j < N_VECTOR; j++) {
            #pragma HLS LOOP_TRIPCOUNT min=N_VECTOR max=N_VECTOR 
            #pragma HLS PIPELINE
            vector_t val = out_stream.read();
            out[(start_i + i) * vec_cols + (vec_start_j + j)] = val;
        }
    }
}


// Matrix-Matrix Multiplication Kernel
extern "C" void krnl_mxm_r(
    data_t *C0, data_t *A0, data_t *B0,
    int start_i, int start_k, int start_j, 
    int M, int N, int P
) {

    vector_t* A = (vector_t*) A0;
    vector_t* B = (vector_t*) B0; 
    vector_t* C = (vector_t*) C0;  

    int vec_start_j = start_j/ VECTOR_SIZE;
    int vec_start_k = start_k/ VECTOR_SIZE;
    int vec_N = N/VECTOR_SIZE;
    int vec_P = P/VECTOR_SIZE;

    #pragma HLS INTERFACE m_axi port=A bundle=gmem0 depth=TILE_SIZE*TILE_SIZE*4
    #pragma HLS INTERFACE m_axi port=B bundle=gmem1 depth=TILE_SIZE*TILE_SIZE*4
    #pragma HLS INTERFACE m_axi port=C bundle=gmem2 depth=TILE_SIZE*TILE_SIZE*4

    #pragma HLS INTERFACE s_axilite port=start_i    bundle=control
    #pragma HLS INTERFACE s_axilite port=start_k    bundle=control
    #pragma HLS INTERFACE s_axilite port=start_j    bundle=control
    #pragma HLS INTERFACE s_axilite port=M          bundle=control
    #pragma HLS INTERFACE s_axilite port=N          bundle=control
    #pragma HLS INTERFACE s_axilite port=P          bundle=control
    #pragma HLS INTERFACE s_axilite port=return     bundle=control

    static stream_t in1("input_stream_A");
    static stream_t in2("input_stream_B");
    static stream_t in3("input_stream_C");

    //static stream_t fromdo_toout("doout_stream");

    static stream_t out("output_stream_C");

    #pragma HLS STREAM variable=in1 depth=32
    #pragma HLS STREAM variable=in2 depth=32
    #pragma HLS STREAM variable=in3 depth=32

    //#pragma HLS STREAM variable=fromdo_toout depth=TILE_SIZE*N_VECTOR

    #pragma HLS STREAM variable=out depth=32

    #pragma HLS dataflow 
        // Load matrix A, B, and C into streams
        load_input_mat(A, in1, start_i, vec_start_k, vec_N);
        load_input_mat(B, in2, start_k, vec_start_j, vec_P);
        load_input_mat(C, in3, start_i, vec_start_j, vec_P);
        //load_input_matricies(A, B, C, in1, in2, in3, start_i, vec_start_j, vec_N, vec_P);

        // Perform matrix-matrix multiplication
        do_calc(out, in1, in2, in3);

        // Store result matrix C
        store_result(C, out, start_i, vec_start_j, vec_P); 

        //accumalte_store_result(C, in3, out, start_i, vec_start_j, vec_P);

}
