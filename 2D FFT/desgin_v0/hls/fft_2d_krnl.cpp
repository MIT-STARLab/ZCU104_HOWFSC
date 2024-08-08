
#include "fft_2d_krnl_include.h"



////////////////////////////////           Data Moving          ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////       

// Data Moving and Control Unit
// void data_mover(cmpx_data_t *input_mat,
//                 cmpx_data_t *output_mat,
//                 cmpx_data_t *temp_mat,
//                 hls::stream<cmpx_data_t> &matrix_to_row_strm,
//                 hls::stream<cmpx_data_t> &matrix_from_row_strm, 
//                 hls::stream<cmpx_data_t> &matrix_to_col_strm,
//                 hls::stream<cmpx_data_t> &matrix_from_col_strm) {

//     // #pragma HLS DATAFLOW
//     stream_to_fft_row(input_mat, matrix_to_row_strm);
//     stream_from_fft_row_to_fft_col(temp_mat, matrix_from_row_strm, matrix_to_col_strm);
//     stream_from_fft_col(output_mat, matrix_from_col_strm);
// }




void stream_to_fft_row(cmpx_data_t *input_mat, hls::stream<cmpx_data_t> &matrix_to_row_strm){

    STRAM_TO_FFT_ROW:
    for (int i=0; i<MAT_SIZE; i++){               // MAT_SIZE = MAT_ROWS * MAT_COLS
        #pragma HLS PIPELINE II=1
        matrix_to_row_strm << input_mat[i];       // matrix is in row major representaion
    }

}


void stream_from_fft_row_to_fft_col(cmpx_data_t *temp_mat,  hls::stream<cmpx_data_t> &matrix_from_row_strm,  hls::stream<cmpx_data_t> &matrix_to_col_strm){

    STRAM_FROM_FFT_ROW:
    for (int i=0; i<MAT_SIZE; i++){
        #pragma HLS PIPELINE II=1
        temp_mat[i] = matrix_from_row_strm.read();
    }

    // The first loop needs to be done before the second loop starts
    // Sequential ORDER MATTER HERE
    // We should only start streaming to col wise fft after finishing the row wise fft

    STRAM_TO_FFT_COL:
    for (int j=0; j<MAT_COLS; j++){
        for (int i=0; i<MAT_ROWS; i++){
            
            #pragma HLS PIPELINE II=1
            matrix_to_col_strm << temp_mat[i*MAT_COLS + j];   // stream the matrix in col major
        }
    }
    
    
}



void stream_from_fft_col(cmpx_data_t *output_mat, hls::stream<cmpx_data_t> &matrix_from_col_strm){

    STRAM_FROM_FFT_COL:
    for (int j=0; j<MAT_COLS; j++){
        for (int i=0; i<MAT_ROWS; i++){
            #pragma HLS PIPELINE II=1
            output_mat[i*MAT_COLS + j] = matrix_from_col_strm.read();
        }
    }
    

}





////////////////////////////////           Row-Wise FFT          ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


void fft_row_top(bool direction, hls::stream<cmpx_data_t> &matrix_to_row_strm, hls::stream<cmpx_data_t> &matrix_from_row_strm){


    LOOP_FFT_ROWS:
    for (int row= 0; row<MAT_ROWS; row++){

        cmpx_data_t row_in[MAT_COLS]__attribute__((no_ctor));
        cmpx_data_t row_out[MAT_COLS]__attribute__((no_ctor));
    
        #pragma HLS STREAM  variable=row_in      depth=MAT_COLS
        #pragma HLS STREAM  variable=row_out     depth=MAT_COLS
        #pragma HLS DATAFLOW


        read_row_in(matrix_to_row_strm, row_in);
        fft_row_unit(direction, row_in, row_out);
        write_row_out(matrix_from_row_strm, row_out);
    }

}



void fft_row_unit(bool direction, cmpx_data_t row_in[MAT_COLS], cmpx_data_t row_out[MAT_COLS]){

    #pragma HLS DATAFLOW
    configRow_t config_row;
    statusRow_t status_row;
    fft_row_init(direction, &config_row);
    hls::fft<apra_2dfft_config_row>(row_in, row_out, &status_row, &config_row);
}



void fft_row_init(bool direction, configRow_t* config_row){
    config_row->setDir(direction);
}



void read_row_in(hls::stream<cmpx_data_t> &matrix_to_row_strm, cmpx_data_t row_in[MAT_COLS]){

    READ_ROW_IN:
    for (int i=0; i< MAT_COLS; i++){
        #pragma HLS PIPELINE II=1
        row_in[i] = matrix_to_row_strm.read();
    }

}


void write_row_out(hls::stream<cmpx_data_t> &matrix_from_row_strm, cmpx_data_t row_out[MAT_COLS]){
    READ_ROW_OUT:
    for (int i=0; i< MAT_COLS; i++){
        #pragma HLS PIPELINE II=1
        matrix_from_row_strm << row_out[i];
    }

}





////////////////////////////////           Col-Wise FFT          ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// Col Wise FFT 
void fft_col_top(bool direction, hls::stream<cmpx_data_t> &matrix_to_col_strm, hls::stream<cmpx_data_t> &matrix_from_col_strm){



    LOOP_FFT_COLS:
    for (int col= 0; col<MAT_COLS; col++){

        cmpx_data_t col_in[MAT_ROWS]__attribute__((no_ctor));
        cmpx_data_t col_out[MAT_ROWS]__attribute__((no_ctor));

        #pragma HLS STREAM  variable=col_in        depth=MAT_ROWS
        #pragma HLS STREAM  variable=col_out       depth=MAT_ROWS
        
        #pragma HLS DATAFLOW

        read_col_in(matrix_to_col_strm, col_in);
        fft_col_unit(direction, col_in, col_out);
        write_col_out(matrix_from_col_strm, col_out);
    }
}


void fft_col_unit(bool direction, cmpx_data_t col_in[MAT_ROWS], cmpx_data_t col_out[MAT_ROWS]){

    #pragma HLS DATAFLOW
    configCol_t config_col;
    statusCol_t status_col;

    fft_col_init(direction, &config_col);
    hls::fft<apra_2dfft_config_col>(col_in, col_out, &status_col, &config_col);

}





void fft_col_init(bool direction, configCol_t* config_col){
    config_col -> setDir(direction);
}


void read_col_in(hls::stream<cmpx_data_t> &matrix_to_col_strm,cmpx_data_t col_in[MAT_ROWS]){

    READ_COL_IN:
    for (int i=0; i< MAT_ROWS; i++){
        #pragma HLS PIPELINE II=1
        col_in[i] =  matrix_to_col_strm.read();
    }


}


void write_col_out(hls::stream<cmpx_data_t> &matrix_from_col_strm, cmpx_data_t col_out[MAT_ROWS]){
    READ_COL_OUT:
    for (int i=0; i< MAT_ROWS; i++){
        #pragma HLS PIPELINE II=1
        matrix_from_col_strm << col_out[i];
    }

}









////////////////////////////////           TOP LEVEL 2D FFT          ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


void fft_2d(bool direction, cmpx_data_t *input_mat, cmpx_data_t *temp_mat, cmpx_data_t *output_mat){

    #pragma HLS INTERFACE s_axilite     port=direction
    #pragma HLS INTERFACE m_axi         port=input_mat      bundle=gmem0  depth = MAT_SIZE
    #pragma HLS INTERFACE m_axi         port=output_mat     bundle=gmem0  depth = MAT_SIZE
    #pragma HLS INTERFACE m_axi         port=temp_mat     bundle=gmem1  depth = MAT_SIZE

    hls::stream<cmpx_data_t> matrix_to_row_strm;
    hls::stream<cmpx_data_t> matrix_from_row_strm; 
    hls::stream<cmpx_data_t> matrix_to_col_strm;
    hls::stream<cmpx_data_t> matrix_from_col_strm;
    // TODO: STREAM PRAGMAS: AXI OR FIFO

    #pragma HLS DATAFLOW
    // data_mover(input_mat, output_mat, temp_mat, matrix_to_row_strm, matrix_from_row_strm, matrix_to_col_strm, matrix_from_col_strm);

    stream_to_fft_row(input_mat, matrix_to_row_strm);
    fft_row_top(direction, matrix_to_row_strm, matrix_from_row_strm);
    stream_from_fft_row_to_fft_col(temp_mat, matrix_from_row_strm,  matrix_to_col_strm);
    fft_col_top(direction, matrix_to_col_strm, matrix_from_col_strm);
    stream_from_fft_col(output_mat, matrix_from_col_strm);
}