/*
 * MIT STAR Lab
 * Last Modified by Subhi in Nov 21, 2024
 * Pure software implementations of useful FFT function to help test hardware kernels
 */

#include "fft.h"

/**
 * Inherts docs
 */
void fft(cmpx_data_t *a, int length, bool invert, bool scale){

    // Fast Failing Assertions
    if (log2(length) != int(log2(length))){
        throw runtime_error("Data array length must be a power of 2");
    }

    int n = length;

    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = - 2 * M_PI / len * (invert ? -1 : 1);
        cmpx_data_t wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            cmpx_data_t w(1);
            for (int j = 0; j < len / 2; j++) {
                cmpx_data_t u = a[i+j], v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }

    // If we want to scale the output of the inverse fft by length
    if (invert && scale) {
        for (int i =0; i< length; i++)
            a[i] /= n;
    }

}

/**
 * Inherts docs
 */
void fft2d(cmpx_data_t *a, int m, int n, bool invert, bool scale){

    // Fast Failing Assertions
    if (log2(n) != int(log2(n))){
        throw runtime_error("n must be a power of 2");
    }

    if (log2(m) != int(log2(m))){
        throw runtime_error("m must be a power of 2");
    }

    //ROWWISE FFT
    for (int row = 0; row < m; row ++){
        cmpx_data_t *row_arr = a + row*n;
        fft(row_arr, n, invert, false);
    }

    //COLWISE FFT
    for (int col = 0; col < n; col++){

        int length = m;
        for (int i = 1, j = 0; i < length; i++) {
            int bit = length >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;

            int index_i = i*n + col;
            int index_j = j*n + col;

            if (i < j)
                swap(a[index_i], a[index_j]);
        }

        for (int len = 2; len <= length; len <<= 1) {
            double ang = - 2 * M_PI / len * (invert ? -1 : 1);
            cmpx_data_t wlen(cos(ang), sin(ang));
            for (int i = 0; i < length; i += len) {
                cmpx_data_t w(1);
                for (int j = 0; j < len / 2; j++) {
                    int i_j_index = (i+j)       *n + col;
                    int i_j_len_2 = (i+j+len/2) *n + col;
                    cmpx_data_t u = a[i_j_index], v = a[i_j_len_2] * w;
                    a[i_j_index] = u + v;
                    a[i_j_len_2] = u - v;
                    w *= wlen;
                }
            }
        }
    }

    // If we want to scale the output of the inverse fft by length
    if (invert && scale) {
        for (int i =0; i< n*m; i++)
            a[i] /= (m*n);
    }

}


/**
 * Inherts docs
 */
void fftshift2d(cmpx_data_t *a, int m, int n)
{
    int half_m = m/2;
    int half_n = n/2;

    int index1, index2, index3, index4;
    int i, j;
    cmpx_data_t tmp;

    for (i = 0; i < half_m; i++) {
        for (j = 0; j < half_n; j++) {

            index1 = i * n + j;
            index3 = (i + half_m) * n + (j+half_n);

            tmp = a[index1];
            a[index1] = a[index3];
            a[index3] = tmp;


            index2 = i * n + (j+half_n);
            index4 = (i + half_m) * n + j;

            tmp = a[index2];
            a[index2] = a[index4];
            a[index4] = tmp;
        }
    }

}


/**
 * Inherts docs
 */
void ifftshift2d(cmpx_data_t *a, int m, int n)
{
    int half_m = m/2;
    int half_n = n/2;

    int index1, index2, index3, index4;
    int i, j;
    cmpx_data_t tmp;

    for (i = 0; i < half_m; i++) {
        for (j = 0; j < half_n; j++) {

            index1 = i * n + j;
            index3 = (i + half_m) * n + (j+half_n);

            tmp = a[index1];
            a[index1] = a[index3];
            a[index3] = tmp;


            index2 = i * n + (j+half_n);
            index4 = (i + half_m) * n + j;

            tmp = a[index2];
            a[index2] = a[index4];
            a[index4] = tmp;
        }
    }

}



//////////////////////////////////    Other Interfaces     //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Inherts docs
 */
void fft_vector(vector<cmpx_data_t> &a, bool invert, bool scale) {

    // Fast Failing Assertions
    if (log2(double(a.size())) != int(log2(double(a.size())))){
        throw runtime_error("Data vector length must be a power of 2");
    }

    int n = a.size();

    // Make reverse-bit order, order of recursion
    // 0 1 2 3 4 5 6 7 -> 0 4 | 2 6 | 1 5 | 3 7
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = - 2 * M_PI / len * (invert ? -1 : 1);
        cmpx_data_t wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            cmpx_data_t w(1);
            for (int j = 0; j < len / 2; j++) {
                cmpx_data_t u = a[i+j], v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }

    // If we want to scale the output of the inverse fft by length
    if (invert && scale) {
        for (cmpx_data_t & x : a)
            x /= n;
    }
}


/**
 * Inherts docs
 */
void fft_recursive(vector<cmpx_data_t> &a, bool invert, bool scale){

    // Fast Failing Assertions
    if (log2(double(a.size())) != int(log2(double(a.size())))){
        throw runtime_error("Data vector length must be a power of 2");
    }

    int n = a.size();

    if (n == 1)
        return;
    
    vector<cmpx_data_t> a0(n/2), a1(n/2);

    for (int i =0; 2*i<n; i++){
        a0[i] = a[2*i];
        a1[i] = a[2*i+1];
    }

    fft_recursive(a0, invert, scale);
    fft_recursive(a1, invert, scale);

    double ang = - 2 * M_PI / n * (invert? -1:1);

    cmpx_data_t w(1), wn(cos(ang), sin(ang));

    for (int i =0; 2*i < n; i++){
        a[i] = a0[i] + w * a1[i];
        a[i + n/2] = a0[i] - w * a1[i];

        // If we want to scale the output of the inverse fft by length
        if (invert && scale){
            a[i]/= 2;
            a[i + n/2] /= 2;
        }

        w *= wn;
    }
}


/**
 * Inherts docs
 */
void fft2d_vector(vector<vector<cmpx_data_t>> &a, bool invert, bool scale) {
    int m = a.size();
    int n = (a[0]).size();

    // Fast Failing Assertions
    if (log2(double(n)) != int(log2(double(n)))){
        throw runtime_error("Data vector length n must be a power of 2");
    }
    if (log2(double(m)) != int(log2(double(m)))){
        throw runtime_error("Data vector length m  must be a power of 2");
    }

    // DO 1D FFT on each row first
    for (int row_index = 0; row_index< m; row_index ++){
        fft_vector(a[row_index], invert, scale);
    }

    // DO 1D FFT on the resulting array
    for (int col_index = 0; col_index < n; col_index ++){
        // cpy col to a vector
        vector<cmpx_data_t> col;
        for (int i = 0; i<m; i++){
            col.push_back(a[i][col_index]);
        }
        // call fft on col
        fft_vector(col, invert, scale);

        // copy back to matrix
        for (int i = 0; i<m; i++){
            a[i][col_index] = col[i];
        }
    }
}




//////////////////////////////    Data Generators    //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


/**
 * Inherts docs
 */
void fft_data_generator_vector(vector<cmpx_data_t> &input_data, vector<cmpx_data_t> &output_data, vector<float> &frequencies, vector<float> &amplitudes, vector<float> &time_range, bool invert, bool scale){

    // Fast Failing Assertions
    if (input_data.size() != output_data.size()){
        throw runtime_error("Input data vector length must match output data vector length");
    }

    if (log2(double(input_data.size())) != int(log2(double(input_data.size())))){
        throw runtime_error("Data vector length must be a power of 2");
    }

    if (frequencies.size() != amplitudes.size()){
        throw runtime_error("frequencies vector size must match vector amplitudes");
    }

    if (time_range.size() != 2){
        throw runtime_error("time_range vector must have only two elements, starting point and an ending point");
    }


    const int fft_data_size = input_data.size();
    const int n = frequencies.size();

    const float time_start = time_range[0];
    const float time_end = time_range[1];
    const float time_interval = (time_end - time_start) / fft_data_size;

    // Fill input data 
    for (int i = 0; i < fft_data_size; i++){
        const float t = time_start + i * time_interval;

        for (int j=0; j<n; j++){
            input_data[i] += (amplitudes[j] * sinf(2 * M_PI * frequencies[j] * t));
            output_data[i] += (amplitudes[j] * sinf(2 * M_PI * frequencies[j] * t));       // we this run the fft function on this as it works inplace
        } 

    }

    // Produce the fft results
    fft_vector(output_data, invert, scale);
}

/**
 * Inherts docs
 */
void fft_data_generator(cmpx_data_t *input_data, cmpx_data_t *output_data, int fft_data_size, vector<float> &frequencies, vector<float> &amplitudes, vector<float> &time_range, bool invert, bool scale){

    if (log2(fft_data_size) != int(log2(fft_data_size))){
        throw runtime_error("Data vector length must be a power of 2");
    }

    if (frequencies.size() != amplitudes.size()){
        throw runtime_error("frequencies vector size must match vector amplitudes");
    }

    if (time_range.size() != 2){
        throw runtime_error("time_range vector must have only two elements, starting point and an ending point");
    }

    const int n = frequencies.size();

    const float time_start = time_range[0];
    const float time_end = time_range[1];
    const float time_interval = (time_end - time_start) / fft_data_size;

    // Fill input data 
    for (int i = 0; i < fft_data_size; i++){
        const float t = time_start + i * time_interval;

        input_data[i] = 0;
        output_data[i] = 0;

        for (int j=0; j<n; j++){
            input_data[i] += (amplitudes[j] * sinf(2 * M_PI * frequencies[j] * t));
            output_data[i] += (amplitudes[j] * sinf(2 * M_PI * frequencies[j] * t));       // we this run the fft function on this as it works inplace
        } 
    }

    // Produce the fft results
    fft(output_data, fft_data_size , invert, scale);
}


/**
 * Inherts docs
 */
void fft_random_data_generator(cmpx_data_t *input_data, cmpx_data_t *output_data, int n, bool invert, bool scale) {

    // Fast Failing Assertions
    if (log2(n) != int(log2(n))){
        throw runtime_error("n must be a power of 2");
    }

    std::random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(100.0, 1000.0);

    // crease random testing data
    for (int i = 0; i < n; ++i) {
        float real = dis(gen);
        float imag = dis(gen);
        input_data [i] = cmpx_data_t(real, imag);
        output_data[i] = cmpx_data_t(real, imag);
    }

    fft(output_data, n, invert, scale);
}





/**
 * Inherts docs
 */
void fft2d_random_data_generator_vector(vector<vector<cmpx_data_t>> &input_data, vector<vector<cmpx_data_t>> &output_data, bool invert, bool scale){
    int m = input_data.size();
    int n = (input_data[0]).size();

    // Fast Failing Assertions
    if (input_data.size() != output_data.size()){
        throw runtime_error("Input data vector length must match output data vector length");
    }

    if (log2(double(input_data.size())) != int(log2(double(input_data.size())))){
        throw runtime_error("Data vector length must be a power of 2");
    }

    std::random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 10.0);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            double real = dis(gen);
            double imag = dis(gen);
            input_data[i][j] = cmpx_data_t(real, imag);
            output_data[i][j] = cmpx_data_t(real, imag);
        }
    }

    fft2d_vector(output_data, invert, scale);
}


/**
 * Inherts docs
 */
void fft2d_random_data_generator(cmpx_data_t *input_data, cmpx_data_t *output_data, int m, int n, bool invert, bool scale) {

    // Fast Failing Assertions
    if (log2(n) != int(log2(n))){
        throw runtime_error("n must be a power of 2");
    }

    if (log2(m) != int(log2(m))){
        throw runtime_error("m must be a power of 2");
    }

    std::random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(100.0, 1000.0);

    // crease random testing data
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int index = i * n + j;
            float real = dis(gen);
            float imag = dis(gen);
            input_data [index] = cmpx_data_t(real, imag);
            output_data[index] = cmpx_data_t(real, imag);
        }
    }

    fft2d(output_data, m, n, invert, scale);
}