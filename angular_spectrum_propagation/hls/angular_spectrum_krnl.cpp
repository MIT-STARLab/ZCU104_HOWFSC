/*
 * MIT STAR Lab
 * Modified by Subhi in Nov 26, 2024
 */

#include "angular_spectrum_include.h"


void shifted_stream_to_first_fft(cmpx_data_t *input_mat, hls::stream<vector_row_data_t> &shifted_matrix_to_first_fft_strm)
{
    vector_row_data_t* input_mat_vec = (vector_row_data_t*) input_mat;

    SHIFTED_STREAM_TO_FIRST_FFT_UPPER_HALF:
        for (int i= HALF_MAT_SIZE/VECTOR_SIZE_ROW; i<MAT_SIZE/VECTOR_SIZE_ROW; i++){


            // half a row from quad 4
            SHIFTED_STREAM_TO_FIRST_FFT_QUAD4:
            for (int j=HALF_MAT_ROWS/VECTOR_SIZE_ROW; j<MAT_ROWS/VECTOR_SIZE_ROW; j++){
                #pragma HLS PIPELINE II=1
                shifted_matrix_to_first_fft_strm << input_mat_vec[i+j];
            }

            SHIFTED_STREAM_TO_FIRST_FFT_QUAD3:
            // then, half a row from quad 3
            for (int j=0; j<HALF_MAT_ROWS/VECTOR_SIZE_ROW; j++){
                #pragma HLS PIPELINE II=1
                shifted_matrix_to_first_fft_strm << input_mat_vec[i+j];
            }
        }

    SHIFTED_STREAM_TO_FIRST_FFT_LOWER_HALF:
        for (int i= 0; i<HALF_MAT_SIZE/VECTOR_SIZE_ROW; i++){


            // half a row from quad 2
            SHIFTED_STREAM_TO_FIRST_FFT_QUAD2:
            for (int j=HALF_MAT_ROWS/VECTOR_SIZE_ROW; j<MAT_ROWS/VECTOR_SIZE_ROW; j++){
                #pragma HLS PIPELINE II=1
                shifted_matrix_to_first_fft_strm << input_mat_vec[i+j];
            }

            SHIFTED_STREAM_TO_FIRST_FFT_QUAD1:
            // then, half a row from quad 1
            for (int j=0; j<HALF_MAT_ROWS/VECTOR_SIZE_ROW; j++){
                #pragma HLS PIPELINE II=1
                shifted_matrix_to_first_fft_strm << input_mat_vec[i+j];
            }
        }
}


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
        #pragma HLS PIPELINE II=1
        matrix_to_row_strm << temp_mat_row_col[i];
    }
}


void shifted_resolve_stream_from_second_fft(cmpx_data_t *output_mat, hls::stream<vector_row_data_t> &second_output_matrix_row_major_strm)
{
    vector_row_data_t* output_mat_vec = (vector_row_data_t*) output_mat;

    SHIFTED_STREAM_BACK_FROM_SECOND_FFT_LOWER_HALF:
    for (int i= HALF_MAT_SIZE/VECTOR_SIZE_ROW; i<MAT_SIZE/VECTOR_SIZE_ROW; i++){

        RESOLVE_SHIFTED_STREAM_FROM_SECOND_FFT_QUAD4:
        for (int j=HALF_MAT_ROWS/VECTOR_SIZE_ROW; j<MAT_ROWS/VECTOR_SIZE_ROW; j++){
            #pragma HLS PIPELINE II=1
            output_mat_vec[i+j] = second_output_matrix_row_major_strm.read();
        }

        RESOLVE_SHIFTED_STREAM_FROM_SECOND_FFT_QUAD3:
        for (int j=0; j<HALF_MAT_ROWS/VECTOR_SIZE_ROW; j++){
            #pragma HLS PIPELINE II=1
            output_mat_vec[i+j] = second_output_matrix_row_major_strm.read();
        }
    }

    SHIFTED_STREAM_BACK_FROM_SECOND_FFT_UPPER_HALF:
    for (int i= 0; i<HALF_MAT_SIZE/VECTOR_SIZE_ROW; i++){

        RESOLVE_SHIFTED_STREAM_FROM_SECOND_FFT_QUAD2:
        for (int j=HALF_MAT_ROWS/VECTOR_SIZE_ROW; j<MAT_ROWS/VECTOR_SIZE_ROW; j++){
            #pragma HLS PIPELINE II=1
            output_mat_vec[i+j] = second_output_matrix_row_major_strm.read();
        }

        RESOLVE_SHIFTED_STREAM_FROM_SECOND_FFT_QUAD1:
        for (int j=0; j<HALF_MAT_ROWS/VECTOR_SIZE_ROW; j++){
            #pragma HLS PIPELINE II=1
            output_mat_vec[i+j] = second_output_matrix_row_major_strm.read();
        }
    }
}


void propagate_wave(
    double scale,
    double distance,
    double k_2,
    double *kxy,
    hls::stream<vector_col_data_t> &first_output_matrix_col_major_strm,
    hls::stream<vector_col_data_t> &second_input_matrix_col_major_strm)
{
    // TODO: LITTLE'S LAW
    for (int col_tile=0; col_tile < HALF_MAT_COLS/VECTOR_SIZE_ROW; col_tile++){
        // this is a quad1 indices that is supposed to be filled with quad4 content
        for (int row=0; row < HALF_MAT_ROWS; row++){
            vector_col_data_t coeff_vec;
            for (int i=0; i<VECTOR_SIZE_COL; i++){
                #pragma HLS UNROLL
                double kxy_i = kxy[row + HALF_MAT_ROWS];
                double kxy_j = kxy[col_tile + i + HALF_MAT_COLS];
                double kxy_sum = kxy_i * kxy_i + kxy_j + kxy_j;
                double kz = hls::sqrt(k_2 - kxy_sum);
                complex<double> e_kzd(1, hls::sin(kz * distance));
                coeff_vec[i] = (cmpx_data_t) (scale * e_kzd);
            }
            second_input_matrix_col_major_strm << coeff_vec * first_output_matrix_col_major_strm.read();
        }

        // this is a quad3 indices that is supposed to be filled with quad2 content
        for (int row=HALF_MAT_ROWS; row < MAT_ROWS; row++){
            vector_col_data_t coeff_vec;
            for (int i=0; i<VECTOR_SIZE_COL; i++){
                #pragma HLS UNROLL
                double kxy_i = kxy[row];
                double kxy_j = kxy[col_tile + i + HALF_MAT_COLS];
                double kxy_sum = kxy_i * kxy_i + kxy_j + kxy_j;
                double kz = hls::sqrt(k_2 - kxy_sum);
                complex<double> e_kzd(1, hls::sin(kz * distance));
                coeff_vec[i] = (cmpx_data_t) (scale * e_kzd);
            }
            second_input_matrix_col_major_strm << coeff_vec * first_output_matrix_col_major_strm.read();
        }
    }


    for (int col_tile=HALF_MAT_COLS/VECTOR_SIZE_ROW; col_tile <MAT_COLS/VECTOR_SIZE_ROW; col_tile++){
        // this is a quad2 indices that is supposed to be filled with quad3 content
        for (int row=0; row < HALF_MAT_ROWS; row++){
            vector_col_data_t coeff_vec;
            for (int i=0; i<VECTOR_SIZE_COL; i++){
                #pragma HLS UNROLL
                double kxy_i = kxy[row + HALF_MAT_ROWS];
                double kxy_j = kxy[col_tile + i];
                double kxy_sum = kxy_i * kxy_i + kxy_j + kxy_j;
                double kz = hls::sqrt(k_2 - kxy_sum);
                complex<double> e_kzd(1, hls::sin(kz * distance));
                coeff_vec[i] = (cmpx_data_t) (scale * e_kzd);
            }
            second_input_matrix_col_major_strm << coeff_vec * first_output_matrix_col_major_strm.read();
        }


        // this is a quad4 indices that is supposed to be filled with quad1 content
        for (int row=HALF_MAT_ROWS; row < MAT_ROWS; row++){
            vector_col_data_t coeff_vec;
            for (int i=0; i<VECTOR_SIZE_COL; i++){
                #pragma HLS UNROLL
                double kxy_i = kxy[row];
                double kxy_j = kxy[col_tile + i];
                double kxy_sum = kxy_i * kxy_i + kxy_j + kxy_j;
                double kz = hls::sqrt(k_2 - kxy_sum);
                complex<double> e_kzd(1, hls::sin(kz * distance));
                coeff_vec[i] = (cmpx_data_t) (scale * e_kzd);
            }
            second_input_matrix_col_major_strm << coeff_vec * first_output_matrix_col_major_strm.read();

        }
    }
}


void angular_spectrum(
    bool direction,
    double scale,
    double distance,
    double k_2,
    double *kxy,
    cmpx_data_t *input_mat,
    cmpx_data_t *output_mat,
    cmpx_data_t *temp_mat_1, 
    cmpx_data_t *temp_mat_2
    )
{
    #pragma HLS INTERFACE s_axilite     port=direction
    #pragma HLS INTERFACE s_axilite     port=scale
    #pragma HLS INTERFACE s_axilite     port=distance
    #pragma HLS INTERFACE s_axilite     port=k_2

    #pragma HLS INTERFACE m_axi         port=kxy              bundle=gmem3  depth = MAT_ROWS        //TODO: cache/BRAM
    #pragma HLS INTERFACE m_axi         port=input_mat        bundle=gmem0  depth = MAT_SIZE
    #pragma HLS INTERFACE m_axi         port=output_mat       bundle=gmem0  depth = MAT_SIZE
    #pragma HLS INTERFACE m_axi         port=temp_mat_1       bundle=gmem1  depth = MAT_SIZE
    #pragma HLS INTERFACE m_axi         port=temp_mat_2       bundle=gmem2  depth = MAT_SIZE

    hls::stream<vector_row_data_t> first_input_matrix_row_major_strm;
    hls::stream<vector_row_data_t> first_output_matrix_row_major_strm;
    hls::stream<vector_col_data_t> first_input_matrix_col_major_strm;
    hls::stream<vector_col_data_t> first_output_matrix_col_major_strm;

    hls::stream<vector_col_data_t> second_input_matrix_col_major_strm;
    hls::stream<vector_col_data_t> second_output_matrix_col_major_strm;
    hls::stream<vector_row_data_t> second_input_matrix_row_major_strm;
    hls::stream<vector_row_data_t> second_output_matrix_row_major_strm;


    #pragma HLS DATAFLOW
    shifted_stream_to_first_fft(input_mat, first_input_matrix_row_major_strm);

    stream_from_fft_row_to_fft_col(temp_mat_1, first_output_matrix_row_major_strm, first_input_matrix_col_major_strm);
    fft_2d(direction, 
           first_input_matrix_row_major_strm,
           first_output_matrix_row_major_strm,
           first_input_matrix_col_major_strm,
           first_output_matrix_col_major_strm);

    propagate_wave(scale, distance, k_2, kxy, first_output_matrix_col_major_strm, second_input_matrix_col_major_strm);

    stream_from_fft_col_to_fft_row(temp_mat_2, second_output_matrix_col_major_strm, second_input_matrix_row_major_strm);
    fft_2d(!direction, //second fft is inverse
           second_input_matrix_row_major_strm,
           second_output_matrix_row_major_strm,
           second_input_matrix_col_major_strm,
           second_output_matrix_col_major_strm);

    shifted_resolve_stream_from_second_fft(output_mat, second_output_matrix_row_major_strm);
}



////////////////////////////////           FFT          /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void fft_2d(
    bool direction, 
    hls::stream<vector_row_data_t> &input_matrix_row_major_strm,
    hls::stream<vector_row_data_t> &output_matrix_row_major_strm,
    hls::stream<vector_col_data_t> &input_matrix_col_major_strm,
    hls::stream<vector_col_data_t> &output_matrix_col_major_strm
    )
{
    #pragma HLS INTERFACE axis port = input_matrix_row_major_strm;
    #pragma HLS INTERFACE axis port = output_matrix_row_major_strm;
    #pragma HLS INTERFACE axis port = input_matrix_col_major_strm;
    #pragma HLS INTERFACE axis port = output_matrix_col_major_strm;

    #pragma HLS DATAFLOW
    fft_row_top(direction, input_matrix_row_major_strm, output_matrix_row_major_strm);
    fft_col_top(direction, input_matrix_col_major_strm, output_matrix_col_major_strm);
}




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


void serialize_rows_in(hls::stream<vector_row_data_t> &matrix_to_row_strm, cmpx_stream_t &row_in1, cmpx_stream_t &row_in2)
{
    SERIALIZE_ROW_IN1:
    for (int i=0; i< MAT_COLS/VECTOR_SIZE_ROW; i++){
        #pragma HLS PIPELINE

        vector_row_data_t row_segment1 = matrix_to_row_strm.read();
        for (int j = 0; j < VECTOR_SIZE_ROW; j++){
            #pragma HLS UNROLL
            row_in1 << row_segment1[j];
        }
        
    }

    SERIALIZE_ROW_IN2:
    for (int i=0; i< MAT_COLS/VECTOR_SIZE_ROW; i++){
        #pragma HLS PIPELINE
        vector_row_data_t row_segment2 = matrix_to_row_strm.read();
        for (int j = 0; j < VECTOR_SIZE_ROW; j++){
            #pragma HLS UNROLL
            row_in2 << row_segment2[j];
        }
        
    }
}

void serialize_rows_out(hls::stream<vector_row_data_t> &matrix_from_row_strm, cmpx_stream_t &row_out1, cmpx_stream_t &row_out2)
{
    SERIALIZE_ROW_OUT1:
    for (int i=0; i< MAT_COLS/VECTOR_SIZE_ROW; i++){

        #pragma HLS PIPELINE II=2
        #pragma pragma HLS DATAFLOW
    
        vector_row_data_t row1_segment;
        for (int j = 0; j < VECTOR_SIZE_ROW; j++){
            // #pragma HLS UNROLL
            row1_segment[j] = row_out1.read();
        }
        matrix_from_row_strm << row1_segment;
    }


    SERIALIZE_ROW_OUT2:
    for (int i=0; i< MAT_COLS/VECTOR_SIZE_ROW; i++){

        #pragma HLS PIPELINE II=2
        #pragma pragma HLS DATAFLOW
    
        vector_row_data_t row2_segment;
        for (int j = 0; j < VECTOR_SIZE_ROW; j++){
            // #pragma HLS UNROLL
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
        #pragma pragma HLS DATAFLOW
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
        #pragma pragma HLS DATAFLOW
        vector_col_data_t col_segment;
        col_segment[0] = col_out1.read();
        col_segment[1] = col_out2.read();
        col_segment[2] = col_out3.read();
        col_segment[3] = col_out4.read();

        matrix_from_col_strm <<  col_segment;
    }
}
