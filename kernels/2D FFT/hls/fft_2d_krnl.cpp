
#include "fft_2d_krnl_include.h"



////////////////////////////////           Data Moving          ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////       


void stream_to_fft_row(cmpx_data_t *input_mat, hls::stream<vector_row_data_t> &matrix_to_row_strm){

    vector_row_data_t* input_mat_vec = (vector_row_data_t*) input_mat;

    STRAM_TO_FFT_ROW:
    for (int i=0; i<MAT_SIZE/VECTOR_SIZE_ROW; i++){               // (MAT_ROWS/ VECTOR_SIZE_ROW) * MAT_COLS
        #pragma HLS PIPELINE II=1
        matrix_to_row_strm << input_mat_vec[i];                       // matrix is in row major representaion
    }

}


void stream_from_fft_row_to_fft_col(cmpx_data_t *temp_mat,  hls::stream<vector_row_data_t> &matrix_from_row_strm,  hls::stream<vector_col_data_t> &matrix_to_col_strm)
{
    vector_row_data_t* temp_map_vec_row = (vector_row_data_t*) temp_mat;

    STRAM_FROM_FFT_ROW:
    for (int i=0; i<MAT_SIZE/VECTOR_SIZE_ROW; i++){

        #pragma HLS PIPELINE II=1
        temp_map_vec_row[i] = matrix_from_row_strm.read();
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



void stream_from_fft_col(cmpx_data_t *output_mat, hls::stream<vector_col_data_t> &matrix_from_col_strm)
{

    vector_col_data_t* output_mat_vec_col = (vector_col_data_t*) output_mat;

    STRAM_FROM_FFT_COL:
    for (int j=0; j<MAT_COLS/VECTOR_SIZE_COL; j++){
        for (int i=0; i<MAT_ROWS; i++){
            #pragma HLS PIPELINE II=1

            output_mat_vec_col[i*(MAT_COLS/VECTOR_SIZE_COL) + j] = matrix_from_col_strm.read();
        }
    }

}





////////////////////////////////           Row-Wise FFT          ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


void fft_row_top(bool direction, hls::stream<vector_row_data_t> &matrix_to_row_strm, hls::stream<vector_row_data_t> &matrix_from_row_strm){

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


void fft_row_unit(bool direction, cmpx_stream_t &row_in, cmpx_stream_t &row_out){

    #pragma HLS DATAFLOW
    hls::stream<configRow_t> config_row;
    hls::stream<statusRow_t> status_row;
    fft_row_init(direction, config_row);
    hls::fft<apra_2dfft_config_row>(row_in, row_out, status_row, config_row);
    status_row.read();
}



void fft_row_init(bool direction, hls::stream<configRow_t> &config_row){
    configRow_t cfg;
    cfg.setDir(direction);
    config_row << cfg;
}


void serialize_rows_in(hls::stream<vector_row_data_t> &matrix_to_row_strm, cmpx_stream_t &row_in1, cmpx_stream_t &row_in2){

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

void serialize_rows_out(hls::stream<vector_row_data_t> &matrix_from_row_strm, cmpx_stream_t &row_out1, cmpx_stream_t &row_out2){

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



// Col Wise FFT 
void fft_col_top(bool direction, hls::stream<vector_col_data_t> &matrix_to_col_strm, hls::stream<vector_col_data_t> &matrix_from_col_strm){

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


void fft_col_unit(bool direction, cmpx_stream_t &col_in, cmpx_stream_t &col_out) {

    #pragma HLS DATAFLOW

    hls::stream<configCol_t> config_col;
    hls::stream<statusCol_t> status_col;

    fft_col_init(direction, config_col);
    hls::fft<apra_2dfft_config_col>(col_in, col_out, status_col, config_col);
    status_col.read();
}



void fft_col_init(bool direction, hls::stream<configCol_t> &config_col){
    configCol_t cfg_c;
    cfg_c.setDir(direction);
    config_col << cfg_c;
}



void distibute_cols_in(hls::stream<vector_col_data_t> &matrix_to_col_strm, cmpx_stream_t &col_in1, cmpx_stream_t &col_in2, cmpx_stream_t &col_in3, cmpx_stream_t &col_in4){
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

void distibute_cols_out(hls::stream<vector_col_data_t> &matrix_from_col_strm, cmpx_stream_t &col_out1, cmpx_stream_t &col_out2, cmpx_stream_t &col_out3, cmpx_stream_t &col_out4){
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




////////////////////////////////           TOP LEVEL 2D FFT          ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


void fft_2d(bool direction, cmpx_data_t *input_mat, cmpx_data_t *temp_mat, cmpx_data_t *output_mat){

    #pragma HLS INTERFACE s_axilite     port=direction
    #pragma HLS INTERFACE m_axi         port=input_mat      bundle=gmem0  depth = MAT_SIZE
    #pragma HLS INTERFACE m_axi         port=output_mat     bundle=gmem0  depth = MAT_SIZE
    #pragma HLS INTERFACE m_axi         port=temp_mat     bundle=gmem1  depth = MAT_SIZE

    hls::stream<vector_row_data_t> matrix_to_row_strm;
    hls::stream<vector_row_data_t> matrix_from_row_strm; 
    hls::stream<vector_col_data_t> matrix_to_col_strm;
    hls::stream<vector_col_data_t> matrix_from_col_strm;
    // TODO: STREAM PRAGMAS: AXI OR FIFO?

    #pragma HLS DATAFLOW

    stream_to_fft_row(input_mat, matrix_to_row_strm);
    fft_row_top(direction, matrix_to_row_strm, matrix_from_row_strm);
    stream_from_fft_row_to_fft_col(temp_mat, matrix_from_row_strm,  matrix_to_col_strm);
    fft_col_top(direction, matrix_to_col_strm, matrix_from_col_strm);
    stream_from_fft_col(output_mat, matrix_from_col_strm);
}