/*
 * MIT STAR Lab
 * M.Subhi Abo Rdan (msubhi_a@mit.edu)
 * Last modified on January 23, 2024
 * Pure software implementations of Angular Spectrum Propagation Method to help test hardware kernels
 */

#include "angular_spectrum_krnl.h"


void shifted_stream_to_first_fft(cmpx_data_t *input_mat, hls::stream<vector_row_data_t> &shifted_matrix_to_first_fft_strm)
{
    vector_row_data_t* input_mat_vec = (vector_row_data_t*) input_mat;

    int start_index = ((MAT_SIZE/2)/VECTOR_SIZE_ROW) + (MAT_ROWS/2)/VECTOR_SIZE_ROW;
    int index = start_index;
    bool flag_add = false;
    int run_step = 0;
    int count = 0;

    SHIFTED_STREAM_TO_FIRST_FFT_UPPER_HALF:
    for (int runs_so_far = 0; runs_so_far<MAT_ROWS; runs_so_far+=run_step){
        #pragma HLS pipeline
        #pragma HLS loop_tripcount  min=((MAT_SIZE/2)/VECTOR_SIZE_ROW)  max = ((MAT_SIZE/2)/VECTOR_SIZE_ROW)
        shifted_matrix_to_first_fft_strm << input_mat_vec[index];
        index++;
        count++;
        run_step = 0;
        if (count == (MAT_ROWS/2)/VECTOR_SIZE_ROW){
            run_step = 1;
            count = 0;
            if (flag_add){
                index += MAT_ROWS/VECTOR_SIZE_ROW;
            }
            else {
                index -= MAT_ROWS/VECTOR_SIZE_ROW;
            }
            flag_add = !flag_add;
        }
    }

    start_index = (MAT_ROWS/2)/VECTOR_SIZE_ROW;
    index = start_index;
    flag_add = false;
    run_step = 0;
    count = 0;

    SHIFTED_STREAM_TO_FIRST_FFT_LOWER_HALF:
    for (int runs_so_far = 0; runs_so_far<MAT_ROWS; runs_so_far+=run_step){
        #pragma HLS pipeline
        #pragma HLS loop_tripcount  min=((MAT_SIZE/2)/VECTOR_SIZE_ROW) max = ((MAT_SIZE/2)/VECTOR_SIZE_ROW)
        shifted_matrix_to_first_fft_strm << input_mat_vec[index];
        index++;
        count++;
        run_step = 0;
        if (count == (MAT_ROWS/2)/VECTOR_SIZE_ROW){
            run_step = 1;
            count = 0;
            if (flag_add){
                index += MAT_ROWS/VECTOR_SIZE_ROW;
            }
            else {
                index -= MAT_ROWS/VECTOR_SIZE_ROW;
            }
            flag_add = !flag_add;
        }
    }

}

// void shifted_stream_to_first_fft(cmpx_data_t *input_mat, hls::stream<vector_row_data_t> &shifted_matrix_to_first_fft_strm)
// {
//     vector_row_data_t* input_mat_vec = (vector_row_data_t*) input_mat;

//     SHIFTED_STREAM_TO_FIRST_FFT_UPPER_HALF:
//         for (int i= ((MAT_SIZE/2)/VECTOR_SIZE_ROW); i<(MAT_SIZE/VECTOR_SIZE_ROW); i+= (MAT_ROWS/VECTOR_SIZE_ROW)){

//             // half a row from quad 4
//             SHIFTED_STREAM_TO_FIRST_FFT_QUAD4:
//             for (int j=(MAT_ROWS/2)/VECTOR_SIZE_ROW; j<MAT_ROWS/VECTOR_SIZE_ROW; j++){
//                 #pragma HLS PIPELINE II=1
//                 shifted_matrix_to_first_fft_strm << input_mat_vec[i+j];
//             }

//             SHIFTED_STREAM_TO_FIRST_FFT_QUAD3:
//             // then, half a row from quad 3
//             for (int j=0; j<(MAT_ROWS/2)/VECTOR_SIZE_ROW; j++){
//                 #pragma HLS PIPELINE II=1
//                 shifted_matrix_to_first_fft_strm << input_mat_vec[i+j];
//             }
//         }

//     SHIFTED_STREAM_TO_FIRST_FFT_LOWER_HALF:
//         for (int i= 0; i<(MAT_SIZE/2)/VECTOR_SIZE_ROW; i+= MAT_ROWS/VECTOR_SIZE_ROW){

//             // half a row from quad 2
//             SHIFTED_STREAM_TO_FIRST_FFT_QUAD2:
//             for (int j=(MAT_ROWS/2)/VECTOR_SIZE_ROW; j<MAT_ROWS/VECTOR_SIZE_ROW; j++){
//                 #pragma HLS PIPELINE II=1
//                 shifted_matrix_to_first_fft_strm << input_mat_vec[i+j];
//             }

//             SHIFTED_STREAM_TO_FIRST_FFT_QUAD1:
//             // then, half a row from quad 1
//             for (int j=0; j<(MAT_ROWS/2)/VECTOR_SIZE_ROW; j++){
//                 #pragma HLS PIPELINE II=1
//                 shifted_matrix_to_first_fft_strm << input_mat_vec[i+j];
//             }
//         }
// }


void stream_from_fft_row_to_fft_col(cmpx_data_t *temp_mat,  hls::stream<vector_row_data_t> &matrix_from_row_strm,  hls::stream<vector_col_data_t> &matrix_to_col_strm)
{
    vector_row_data_t* temp_mat_vec_row = (vector_row_data_t*) temp_mat;
    STRAM_FROM_FFT_ROW:
    for (int i=0; i<MAT_SIZE/VECTOR_SIZE_ROW; i++){
        #pragma HLS PIPELINE II=1
        temp_mat_vec_row[i] = matrix_from_row_strm.read();
    }
    vector_col_data_t* temp_mat_vec_col = (vector_col_data_t*) temp_mat;

    STRAM_TO_FFT_COL:
    for (int j=0; j<MAT_COLS/VECTOR_SIZE_COL; j++){
        for (int i=0; i<MAT_ROWS; i++){
            
            #pragma HLS PIPELINE II=1
            matrix_to_col_strm << temp_mat_vec_col[i*(MAT_COLS/VECTOR_SIZE_COL) + j];   // stream the matrix in col major
        }
    }
}


void stream_from_fft_col_to_fft_row(cmpx_data_t *temp_mat,  hls::stream<vector_col_data_t> &matrix_from_col_strm,  hls::stream<vector_row_data_t> &matrix_to_row_strm)
{
    vector_col_data_t* temp_mat_vec_col = (vector_col_data_t*) temp_mat;

    STRAM_FROM_FFT_COL:
    for (int j=0; j<MAT_COLS/VECTOR_SIZE_COL; j++){
        for (int i=0; i<MAT_ROWS; i++){
            #pragma HLS PIPELINE II=1
            temp_mat_vec_col[i*(MAT_COLS/VECTOR_SIZE_COL) + j] = matrix_from_col_strm.read();
        }
    }

    vector_row_data_t* temp_mat_row_col = (vector_row_data_t*) temp_mat;

    STRAM_TO_FFT_ROW:
    for (int i=0; i<MAT_SIZE/VECTOR_SIZE_ROW; i++){               // (MAT_ROWS/ VECTOR_SIZE_ROW) * MAT_COLS
        #pragma HLS PIPELINE
        matrix_to_row_strm << temp_mat_row_col[i];
    }
}


// void shifted_resolve_stream_from_second_fft(cmpx_data_t *output_mat, hls::stream<vector_row_data_t> &second_output_matrix_row_major_strm)
// {
//     vector_row_data_t* output_mat_vec = (vector_row_data_t*) output_mat;

//     SHIFTED_STREAM_BACK_FROM_SECOND_FFT_LOWER_HALF:
//     for (int i= (MAT_SIZE/2)/VECTOR_SIZE_ROW; i<MAT_SIZE/VECTOR_SIZE_ROW; i+= MAT_ROWS/VECTOR_SIZE_ROW){

//         RESOLVE_SHIFTED_STREAM_FROM_SECOND_FFT_QUAD4:
//         for (int j=(MAT_ROWS/2)/VECTOR_SIZE_ROW; j<MAT_ROWS/VECTOR_SIZE_ROW; j++){
//             #pragma HLS PIPELINE II=1
//             output_mat_vec[i+j] = second_output_matrix_row_major_strm.read();
//         }

//         RESOLVE_SHIFTED_STREAM_FROM_SECOND_FFT_QUAD3:
//         for (int j=0; j<(MAT_ROWS/2)/VECTOR_SIZE_ROW; j++){
//             #pragma HLS PIPELINE II=1
//             output_mat_vec[i+j] = second_output_matrix_row_major_strm.read();
//         }
//     }

//     SHIFTED_STREAM_BACK_FROM_SECOND_FFT_UPPER_HALF:
//     for (int i= 0; i<(MAT_SIZE/2)/VECTOR_SIZE_ROW; i+= MAT_ROWS/VECTOR_SIZE_ROW){

//         RESOLVE_SHIFTED_STREAM_FROM_SECOND_FFT_QUAD2:
//         for (int j=(MAT_ROWS/2)/VECTOR_SIZE_ROW; j<MAT_ROWS/VECTOR_SIZE_ROW; j++){
//             #pragma HLS PIPELINE II=1
//             output_mat_vec[i+j] = second_output_matrix_row_major_strm.read();
//         }

//         RESOLVE_SHIFTED_STREAM_FROM_SECOND_FFT_QUAD1:
//         for (int j=0; j<(MAT_ROWS/2)/VECTOR_SIZE_ROW; j++){
//             #pragma HLS PIPELINE II=1
//             output_mat_vec[i+j] = second_output_matrix_row_major_strm.read();
//         }
//     }
// }


void shifted_resolve_stream_from_second_fft(cmpx_data_t *output_mat, hls::stream<vector_row_data_t> &second_output_matrix_row_major_strm)
{
    vector_row_data_t* output_mat_vec = (vector_row_data_t*) output_mat;

    int start_index = ((MAT_SIZE/2)/VECTOR_SIZE_ROW) + (MAT_ROWS/2)/VECTOR_SIZE_ROW;
    int index = start_index;
    bool flag_add = false;
    int run_step = 0;
    int count = 0;

    SHIFTED_STREAM_BACK_FROM_SECOND_FFT_LOWER_HALF:
    for (int runs_so_far = 0; runs_so_far<MAT_ROWS; runs_so_far+=run_step){
        #pragma HLS pipeline
        #pragma HLS loop_tripcount  min=((MAT_SIZE/2)/VECTOR_SIZE_ROW)  max = ((MAT_SIZE/2)/VECTOR_SIZE_ROW)
        
        output_mat_vec[index] = second_output_matrix_row_major_strm.read();

        index++;
        count++;
        run_step = 0;
        if (count == (MAT_ROWS/2)/VECTOR_SIZE_ROW){
            run_step = 1;
            count = 0;
            if (flag_add){
                index += MAT_ROWS/VECTOR_SIZE_ROW;
            }
            else {
                index -= MAT_ROWS/VECTOR_SIZE_ROW;
            }
            flag_add = !flag_add;
        }
    }

    start_index = (MAT_ROWS/2)/VECTOR_SIZE_ROW;
    index = start_index;
    flag_add = false;
    run_step = 0;
    count = 0;

    SHIFTED_STREAM_TO_FIRST_FFT_LOWER_HALF:
    for (int runs_so_far = 0; runs_so_far<MAT_ROWS; runs_so_far+=run_step){
        #pragma HLS pipeline
        #pragma HLS loop_tripcount  min=((MAT_SIZE/2)/VECTOR_SIZE_ROW) max = ((MAT_SIZE/2)/VECTOR_SIZE_ROW)

        output_mat_vec[index] = second_output_matrix_row_major_strm.read();

        index++;
        count++;
        run_step = 0;
        if (count == (MAT_ROWS/2)/VECTOR_SIZE_ROW){
            run_step = 1;
            count = 0;
            if (flag_add){
                index += MAT_ROWS/VECTOR_SIZE_ROW;
            }
            else {
                index -= MAT_ROWS/VECTOR_SIZE_ROW;
            }
            flag_add = !flag_add;
        }
    }
}



/**
 * @note The implementation takes the mod of the exponential argumnet by 2 pi. Because high accuracy 
 * needed the order of magnitude of the integer part is high.
 */
cmpx_data_t compute_tf_phase_element(data_t scale, data_t kxy_i, data_t kxy_j, data_t k_2, data_t distance){
    #pragma HLS INLINE off
    #pragma HLS FUNCTION_INSTANTIATE variable=scale, k_2, distance 
    // #pragma HLS DATAFLOW
    data_t kxy_sum;
    data_t kz2;
    data_t kz;
    data_t kzd;
    data_t kzd_mod_2PI;
    data_t sin_kzd;
    data_t cos_kzd;
    hls::complex<data_t> phase_element;

    // primitivedsp has way less latency but not supported on zcu104
    #pragma HLS BIND_OP  variable=kxy_sum 		    op=fmul 	impl=fulldsp  
    #pragma HLS BIND_OP  variable=kxy_sum 		    op=fadd 	impl=fulldsp  
    #pragma HLS BIND_OP  variable=kz2 			    op=fsub 	impl=fulldsp  
    #pragma HLS BIND_OP  variable=kzd 			    op=fmul 	impl=fulldsp  
    #pragma HLS BIND_OP  variable=phase_element 	op=fmul 	impl=fulldsp  

    const data_t M_2PI = (M_PI * 2);

    kxy_sum = kxy_i * kxy_i + kxy_j * kxy_j;
    kz2 = k_2 - kxy_sum;
    kz = hls::sqrt(kz2);
    kzd = kz * distance;
    kzd_mod_2PI = hls::fmod(kzd, M_2PI);

    sin_kzd = hls::sin(kzd_mod_2PI);
    cos_kzd = hls::cos(kzd_mod_2PI);
    hls::complex<data_t> e_kzd(cos_kzd, sin_kzd);
    phase_element = scale * e_kzd;
    return phase_element;
}

void propagate_wave(
    data_t scale,
    data_t distance,
    data_t k_2,
    data_t *kxy,
    hls::stream<vector_col_data_t> &first_output_matrix_col_major_strm,
    hls::stream<vector_col_data_t> &second_input_matrix_col_major_strm)
{
    // TODO: LITTLE'S LAW; PIPLINING?
    hls::vector<data_t, VECTOR_SIZE_COL> *kxy_col_vectors = (hls::vector<data_t, VECTOR_SIZE_COL> *) kxy;

    int col_tile_index1 = MAT_COLS/(2*VECTOR_SIZE_COL);
    for (int col_tile=MAT_COLS/2; col_tile < MAT_COLS; col_tile+=VECTOR_SIZE_COL){
        hls::vector<data_t, VECTOR_SIZE_COL> kxy_j_vec = kxy_col_vectors[col_tile_index1];

        // this is a quad1 indices that is supposed to be filled with quad4 content
        UPPERLEFT_BUT_QUAD4:
        for (int row=MAT_ROWS/2; row < MAT_ROWS; row++){
            #pragma HLS PIPELINE
            vector_col_data_t coeff_vec;
            data_t kxy_i = kxy[row];

            for (int ee=0; ee<VECTOR_SIZE_COL; ee++){
                #pragma HLS UNROLL
                coeff_vec[ee] = compute_tf_phase_element( scale, kxy_i, kxy_j_vec[ee], k_2, distance);
            }
            second_input_matrix_col_major_strm << coeff_vec * first_output_matrix_col_major_strm.read();
        }

        // this is a quad3 indices that is supposed to be filled with quad2 content
        LOWERLEFT_BUT_QUAD2:
        for (int row=0; row < MAT_ROWS/2; row++){
            #pragma HLS PIPELINE
            vector_col_data_t coeff_vec;
            data_t kxy_i = kxy[row];

            for (int ee=0; ee<VECTOR_SIZE_COL; ee++){
                #pragma HLS UNROLL
                coeff_vec[ee] = compute_tf_phase_element( scale, kxy_i, kxy_j_vec[ee], k_2, distance);
            }
            second_input_matrix_col_major_strm << coeff_vec * first_output_matrix_col_major_strm.read();
        }

        col_tile_index1++;
    }

    int col_tile_index2 = 0;
    for (int col_tile=0; col_tile <MAT_COLS/2; col_tile+=VECTOR_SIZE_COL){
        hls::vector<data_t, VECTOR_SIZE_COL> kxy_j_vec = kxy_col_vectors[col_tile_index2];

        // this is a quad2 indices that is supposed to be filled with quad3 content
        UPPERRIGHT_BUT_QUAD3:
        for (int row=MAT_ROWS/2; row < MAT_ROWS; row++){
            #pragma HLS PIPELINE
            vector_col_data_t coeff_vec;
            data_t kxy_i = kxy[row];

            for (int ee=0; ee<VECTOR_SIZE_COL; ee++){
                #pragma HLS UNROLL
                coeff_vec[ee] = compute_tf_phase_element( scale, kxy_i, kxy_j_vec[ee], k_2, distance);
            }
            second_input_matrix_col_major_strm << coeff_vec * first_output_matrix_col_major_strm.read();
        }

        // this is a quad4 indices that is supposed to be filled with quad1 content
        lowerRIGHT_BUT_QUAD1:
        for (int row=0; row < MAT_ROWS/2; row++){
            #pragma HLS PIPELINE

            vector_col_data_t coeff_vec;
            data_t kxy_i = kxy[row];

            for (int ee=0; ee<VECTOR_SIZE_COL; ee++){
                #pragma HLS UNROLL
                coeff_vec[ee] = compute_tf_phase_element( scale, kxy_i, kxy_j_vec[ee], k_2, distance);
            }
            second_input_matrix_col_major_strm << coeff_vec * first_output_matrix_col_major_strm.read();
        }
        col_tile_index2++;
    }
}


void stream_from_fft_row(cmpx_data_t *output_mat, hls::stream<vector_row_data_t> &first_output_matrix_row_major_strm)
{
    vector_row_data_t* temp_mat_vec_row = (vector_row_data_t*) output_mat;
    STRAM_FROM_FFT_ROW:
    for (int i=0; i<MAT_SIZE/VECTOR_SIZE_ROW; i++){
        #pragma HLS PIPELINE II=16
        temp_mat_vec_row[i] = first_output_matrix_row_major_strm.read();
    }
}

void f1(bool direction, cmpx_data_t *input_mat, cmpx_data_t *output_mat)
{
    hls::stream<vector_row_data_t> first_input_matrix_row_major_strm;
    hls::stream<vector_row_data_t> first_output_matrix_row_major_strm;
    #pragma HLS STREAM  variable=first_input_matrix_row_major_strm    type=fifo
    #pragma HLS STREAM  variable=first_output_matrix_row_major_strm   type=fifo
    #pragma HLS DATAFLOW
    shifted_stream_to_first_fft(input_mat, first_input_matrix_row_major_strm);
    fft_row_top(direction, first_input_matrix_row_major_strm, first_output_matrix_row_major_strm);
    stream_from_fft_row(output_mat, first_output_matrix_row_major_strm);
}


void stream_to_fft_col(cmpx_data_t *output_mat, hls::stream<vector_col_data_t> &first_input_matrix_col_major_strm)
{
    vector_col_data_t* temp_mat_vec_col = (vector_col_data_t*) output_mat;
    STRAM_TO_FFT_COL:
    for (int j=0; j<MAT_COLS/VECTOR_SIZE_COL; j++){
        for (int i=0; i<MAT_ROWS; i++){
            #pragma HLS PIPELINE II=1
            first_input_matrix_col_major_strm << temp_mat_vec_col[i*(MAT_COLS/VECTOR_SIZE_COL) + j];   // stream the matrix in col major
        }
    }
    // int start_index = 0;
    // int index = start_index;
    // int run_step = 0;
    // int count = 0;
    // STRAM_TO_FFT_COL:
    // for (int runs_so_far = 0; runs_so_far<MAT_COLS/VECTOR_SIZE_COL; runs_so_far+=run_step){
    //     #pragma HLS PIPELINE
    //     #pragma HLS loop_tripcount  min=(MAT_SIZE/VECTOR_SIZE_COL) max = (MAT_SIZE/VECTOR_SIZE_COL)
    //     first_input_matrix_col_major_strm << temp_mat_vec_col[index];
    //     index+= (MAT_ROWS/VECTOR_SIZE_COL);
    //     count++;
    //     run_step = 0;
    //     if (count == MAT_ROWS){
    //         run_step = 1;
    //         count = 0;
    //         index -= ((MAT_SIZE/VECTOR_SIZE_ROW) -1);
    //     }
    // }
}



void stream_from_fft_col(cmpx_data_t *input_mat, hls::stream<vector_col_data_t> &second_output_matrix_col_major_strm)
{
    vector_col_data_t* temp_mat_vec_col2 = (vector_col_data_t*) input_mat;

    STRAM_FROM_FFT_COL:
    for (int j=0; j<MAT_COLS/VECTOR_SIZE_COL; j++){
        for (int i=0; i<MAT_ROWS; i++){
            #pragma HLS PIPELINE II=1
            temp_mat_vec_col2[i*(MAT_COLS/VECTOR_SIZE_COL) + j] = second_output_matrix_col_major_strm.read();
        }
    }
    // start_index = 0;
    // index = start_index;
    // run_step = 0;
    // count = 0;

    // STRAM_FROM_FFT_COL:
    // for (int runs_so_far = 0; runs_so_far<MAT_COLS/VECTOR_SIZE_COL; runs_so_far+=run_step){
    //     #pragma HLS PIPELINE
    //     #pragma HLS loop_tripcount  min=(MAT_SIZE/VECTOR_SIZE_COL) max = (MAT_SIZE/VECTOR_SIZE_COL)
    //     temp_mat_vec_col2[index] = second_output_matrix_col_major_strm.read();
    //     index+= (MAT_ROWS/VECTOR_SIZE_COL);
    //     count++;
    //     run_step = 0;
    //     if (count == MAT_ROWS){
    //         run_step = 1;
    //         count = 0;
    //         index -= ((MAT_SIZE/VECTOR_SIZE_ROW) -1);
    //     }
    // }
}

void f2(
        bool direction,
        cmpx_data_t *input_mat,
        cmpx_data_t *output_mat, 
        data_t scale,
        data_t distance,
        data_t k_2,
        data_t *kxy)
{
    bool inv_direction = !direction;
    hls::stream<vector_col_data_t> first_input_matrix_col_major_strm;
    hls::stream<vector_col_data_t> first_output_matrix_col_major_strm;
    hls::stream<vector_col_data_t> second_input_matrix_col_major_strm;
    hls::stream<vector_col_data_t> second_output_matrix_col_major_strm;
    #pragma HLS STREAM  variable=first_input_matrix_col_major_strm    type=fifo
    #pragma HLS STREAM  variable=first_output_matrix_col_major_strm   type=fifo
    #pragma HLS STREAM  variable=second_input_matrix_col_major_strm    type=fifo
    #pragma HLS STREAM  variable=second_output_matrix_col_major_strm   type=fifo
    #pragma HLS DATAFLOW

    stream_to_fft_col(output_mat, first_input_matrix_col_major_strm);
    fft_col_top(direction, first_input_matrix_col_major_strm, first_output_matrix_col_major_strm);
    propagate_wave(scale, distance, k_2, kxy, first_output_matrix_col_major_strm, second_input_matrix_col_major_strm);
    fft_col_top(inv_direction, second_input_matrix_col_major_strm, second_output_matrix_col_major_strm);
    stream_from_fft_col(input_mat, second_output_matrix_col_major_strm);
}


void stream_to_fft_row(cmpx_data_t *input_mat, hls::stream<vector_row_data_t> &second_input_matrix_row_major_strm)
{
    vector_row_data_t* temp_mat_row_col = (vector_row_data_t*) input_mat;
    STRAM_TO_FFT_ROW:
    for (int i=0; i<MAT_SIZE/VECTOR_SIZE_ROW; i++){               // (MAT_ROWS/ VECTOR_SIZE_ROW) * MAT_COLS
        #pragma HLS PIPELINE II=16
        second_input_matrix_row_major_strm << temp_mat_row_col[i];
    }
}

void f3(bool direction, cmpx_data_t *input_mat, cmpx_data_t *output_mat)
{
    bool inv_direction = !direction;
    hls::stream<vector_row_data_t> second_input_matrix_row_major_strm;
    hls::stream<vector_row_data_t> second_output_matrix_row_major_strm;
    
    #pragma HLS STREAM  variable=second_input_matrix_row_major_strm     type=fifo
    #pragma HLS STREAM  variable=second_output_matrix_row_major_strm    type=fifo
    #pragma HLS DATAFLOW
    
    stream_to_fft_row(input_mat, second_input_matrix_row_major_strm);
    fft_row_top(inv_direction, second_input_matrix_row_major_strm, second_output_matrix_row_major_strm);
    shifted_resolve_stream_from_second_fft(output_mat, second_output_matrix_row_major_strm);
}


void angular_spectrum(
    bool direction,
    data_t scale,
    data_t distance,
    data_t k_2,
    data_t *kxy,
    cmpx_data_t *input_mat,
    cmpx_data_t *output_mat
    )
{
    #pragma HLS INTERFACE s_axilite     port=direction
    #pragma HLS INTERFACE s_axilite     port=scale
    #pragma HLS INTERFACE s_axilite     port=distance
    #pragma HLS INTERFACE s_axilite     port=k_2

    #pragma HLS INTERFACE m_axi         port=kxy              bundle=gmem2  depth = MAT_ROWS
    #pragma HLS INTERFACE m_axi         port=input_mat        bundle=gmem0  depth = MAT_SIZE
    #pragma HLS INTERFACE m_axi         port=output_mat       bundle=gmem1  depth = MAT_SIZE

    data_t kxy_uram[MAT_ROWS]; //total size
    #pragma HLS BIND_STORAGE variable=kxy_uram  type=RAM_2P   impl=uram
    for (int i=0; i<MAT_ROWS; i++){ kxy_uram[i] = kxy[i]; }

    #pragma HLS allocation  function    instances=fft_row_top   limit=1
    // Sequential, not DATAFLOW 
    f1(direction,input_mat, output_mat);
    f2(direction,input_mat, output_mat,  scale, distance, k_2, kxy_uram);
    f3(direction,input_mat, output_mat);
}









////////////////////////////////           FFT          /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void fft_row_top(bool direction, hls::stream<vector_row_data_t> &matrix_to_row_strm, hls::stream<vector_row_data_t> &matrix_from_row_strm)
{
    cmpx_stream_t row_in1;
    cmpx_stream_t row_in2;

    cmpx_stream_t row_out1;
    cmpx_stream_t row_out2;

    #pragma HLS STREAM  variable=row_in1      depth=MAT_COLS
    #pragma HLS STREAM  variable=row_in2      depth=MAT_COLS
    #pragma HLS STREAM  variable=row_out1     depth=MAT_COLS
    #pragma HLS STREAM  variable=row_out2     depth=MAT_COLS

    LOOP_FFT_ROWS:
    for (int row= 0; row<MAT_ROWS/FFT_CORE_ROW_COUNT; row++){
        #pragma HLS DATAFLOW
        serialize_rows_in(matrix_to_row_strm, row_in1, row_in2);
        fft_row_unit(direction, row_in1, row_out1);
        fft_row_unit(direction, row_in2, row_out2);
        serialize_rows_out(matrix_from_row_strm, row_out1, row_out2);
    }

}


void fft_row_unit(bool direction, cmpx_stream_t &row_in, cmpx_stream_t &row_out)
{

    #pragma HLS DATAFLOW
    hls::stream<configRow_t> config_row;
    hls::stream<statusRow_t> status_row;
    fft_row_init(direction, config_row);
    hls::fft<apra_2dfft_config_row>(row_in, row_out, status_row, config_row);
    status_row.read();
}



void fft_row_init(bool direction, hls::stream<configRow_t> &config_row)
{
    configRow_t cfg;
    cfg.setDir(direction);
    config_row << cfg;
}


void serialize_rows_in(hls::stream<vector_row_data_t> &matrix_to_row_strm, cmpx_stream_t &row_in1, cmpx_stream_t &row_in2){

    SERIALIZE_ROW_IN1:
    for (int i=0; i< MAT_COLS/VECTOR_SIZE_ROW; i++){
        #pragma HLS PIPELINE II=9

        vector_row_data_t row_segment1 = matrix_to_row_strm.read();
        for (int j = 0; j < VECTOR_SIZE_ROW; j++){
            #pragma HLS UNROLL
            row_in1 << row_segment1[j];
        }
    }

    SERIALIZE_ROW_IN2:
    for (int i=0; i< MAT_COLS/VECTOR_SIZE_ROW; i++){
        #pragma HLS PIPELINE II=9
        vector_row_data_t row_segment2 = matrix_to_row_strm.read();
        for (int j = 0; j < VECTOR_SIZE_ROW; j++){
            #pragma HLS UNROLL
            row_in2 << row_segment2[j];
        }
        
    }
}

void serialize_rows_out(hls::stream<vector_row_data_t> &matrix_from_row_strm, cmpx_stream_t &row_out1, cmpx_stream_t &row_out2){

    SERIALIZE_ROW_OUT1:
    for (int i=0; i< MAT_COLS/VECTOR_SIZE_ROW; i++){

        #pragma HLS PIPELINE II=9
    
        vector_row_data_t row1_segment;
        for (int j = 0; j < VECTOR_SIZE_ROW; j++){
            #pragma HLS UNROLL
            row1_segment[j] = row_out1.read();
        }
        matrix_from_row_strm << row1_segment;
    }


    SERIALIZE_ROW_OUT2:
    for (int i=0; i< MAT_COLS/VECTOR_SIZE_ROW; i++){

        #pragma HLS PIPELINE II=9
    
        vector_row_data_t row2_segment;
        for (int j = 0; j < VECTOR_SIZE_ROW; j++){
            #pragma HLS UNROLL
            row2_segment[j] = row_out2.read();
        }
        matrix_from_row_strm << row2_segment;
    }

}






////////////////////////////////           Col-Wise FFT          ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
void fft_col_top(bool direction, hls::stream<vector_col_data_t> &matrix_to_col_strm, hls::stream<vector_col_data_t> &matrix_from_col_strm)
{
    cmpx_stream_t col_in1;
    cmpx_stream_t col_in2;
    cmpx_stream_t col_in3;
    cmpx_stream_t col_in4;
    cmpx_stream_t col_out1;
    cmpx_stream_t col_out2;
    cmpx_stream_t col_out3;
    cmpx_stream_t col_out4;

    #pragma HLS STREAM  variable=col_in1        depth=MAT_ROWS
    #pragma HLS STREAM  variable=col_in2        depth=MAT_ROWS
    #pragma HLS STREAM  variable=col_in3        depth=MAT_ROWS
    #pragma HLS STREAM  variable=col_in4        depth=MAT_ROWS
    #pragma HLS STREAM  variable=col_out1       depth=MAT_ROWS
    #pragma HLS STREAM  variable=col_out2       depth=MAT_ROWS
    #pragma HLS STREAM  variable=col_out3       depth=MAT_ROWS
    #pragma HLS STREAM  variable=col_out4       depth=MAT_ROWS

    LOOP_FFT_COLS:
    for (int col= 0; col<MAT_COLS/FFT_CORE_COL_COUNT; col++){

        #pragma HLS DATAFLOW
        distibute_cols_in(matrix_to_col_strm, col_in1, col_in2, col_in3, col_in4);

        fft_col_unit(direction, col_in1, col_out1);
        fft_col_unit(direction, col_in2, col_out2);
        fft_col_unit(direction, col_in3, col_out3);
        fft_col_unit(direction, col_in4, col_out4);

        distibute_cols_out(matrix_from_col_strm, col_out1, col_out2, col_out3, col_out4);
    }
}


void fft_col_unit(bool direction, cmpx_stream_t &col_in, cmpx_stream_t &col_out)
{

    #pragma HLS DATAFLOW

    hls::stream<configCol_t> config_col;
    hls::stream<statusCol_t> status_col;

    fft_col_init(direction, config_col);
    hls::fft<apra_2dfft_config_col>(col_in, col_out, status_col, config_col);
    status_col.read();
}



void fft_col_init(bool direction, hls::stream<configCol_t> &config_col)
{
    configCol_t cfg_c;
    cfg_c.setDir(direction);
    config_col << cfg_c;
}



void distibute_cols_in(hls::stream<vector_col_data_t> &matrix_to_col_strm, cmpx_stream_t &col_in1, cmpx_stream_t &col_in2, cmpx_stream_t &col_in3, cmpx_stream_t &col_in4)
{
    DISTRIBUTE_COLS_IN:
    for (int i=0; i< MAT_ROWS; i++){
        #pragma HLS PIPELINE II=1
        vector_col_data_t col_segment = matrix_to_col_strm.read();
        col_in1 << col_segment[0];
        col_in2 << col_segment[1];
        col_in3 << col_segment[2];
        col_in4 << col_segment[3];
    }
}

void distibute_cols_out(hls::stream<vector_col_data_t> &matrix_from_col_strm, cmpx_stream_t &col_out1, cmpx_stream_t &col_out2, cmpx_stream_t &col_out3, cmpx_stream_t &col_out4)
{
    DISTRIBUTE_COLS_OUT:
    for (int i=0; i< MAT_ROWS; i++){
        #pragma HLS PIPELINE II=1
        vector_col_data_t col_segment;
        col_segment[0] = col_out1.read();
        col_segment[1] = col_out2.read();
        col_segment[2] = col_out3.read();
        col_segment[3] = col_out4.read();

        matrix_from_col_strm <<  col_segment;
    }
}





 
/////// Old Design  //////////

// void angular_spectrum(
//     bool direction,
//     data_t scale,
//     data_t distance,
//     data_t k_2,
//     data_t *kxy,
//     cmpx_data_t *input_mat,
//     cmpx_data_t *output_mat,
//     cmpx_data_t *temp_mat_1, 
//     cmpx_data_t *temp_mat_2
//     )
// {
//     #pragma HLS INTERFACE s_axilite     port=direction
//     #pragma HLS INTERFACE s_axilite     port=scale
//     #pragma HLS INTERFACE s_axilite     port=distance
//     #pragma HLS INTERFACE s_axilite     port=k_2

//     //TODO: cache/BRAM
//     #pragma HLS INTERFACE m_axi         port=kxy              bundle=gmem3  depth = MAT_ROWS
//     #pragma HLS INTERFACE m_axi         port=input_mat        bundle=gmem0  depth = MAT_SIZE
//     #pragma HLS INTERFACE m_axi         port=output_mat       bundle=gmem0  depth = MAT_SIZE
//     #pragma HLS INTERFACE m_axi         port=temp_mat_1       bundle=gmem1  depth = MAT_SIZE
//     #pragma HLS INTERFACE m_axi         port=temp_mat_2       bundle=gmem2  depth = MAT_SIZE

//     hls::stream<vector_row_data_t> first_input_matrix_row_major_strm;
//     hls::stream<vector_row_data_t> first_output_matrix_row_major_strm;
//     hls::stream<vector_col_data_t> first_input_matrix_col_major_strm;
//     hls::stream<vector_col_data_t> first_output_matrix_col_major_strm;

//     hls::stream<vector_col_data_t> second_input_matrix_col_major_strm;
//     hls::stream<vector_col_data_t> second_output_matrix_col_major_strm;
//     hls::stream<vector_row_data_t> second_input_matrix_row_major_strm;
//     hls::stream<vector_row_data_t> second_output_matrix_row_major_strm;

//     #pragma HLS DATAFLOW

//     shifted_stream_to_first_fft(input_mat, first_input_matrix_row_major_strm);
//     fft_row_top(direction, first_input_matrix_row_major_strm, first_output_matrix_row_major_strm);
//     stream_from_fft_row_to_fft_col(temp_mat_1, first_output_matrix_row_major_strm, first_input_matrix_col_major_strm);
//     fft_col_top(direction, first_input_matrix_col_major_strm, first_output_matrix_col_major_strm);

//     propagate_wave(scale, distance, k_2, kxy, first_output_matrix_col_major_strm, second_input_matrix_col_major_strm);

//     fft_col_top(!direction, second_input_matrix_col_major_strm, second_output_matrix_col_major_strm);
//     stream_from_fft_col_to_fft_row(temp_mat_2, second_output_matrix_col_major_strm, second_input_matrix_row_major_strm);
//     fft_row_top(!direction, second_input_matrix_row_major_strm, second_output_matrix_row_major_strm);

//     shifted_resolve_stream_from_second_fft(output_mat, second_output_matrix_row_major_strm);
// }
