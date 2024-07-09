
#include "fft.hxx"

/**
 * 1D FFT Implementation taken from https://cp-algorithms.com/algebra/fft.html (accessed on Jun 25).
 */


/**
 * Inherts docs
 */
void fft(vector<cmpx_data_t> &a, bool invert) {

    // Fast Failing Assertions
    if (log2(double(a.size())) != int(log2(double(a.size())))){
        throw runtime_error("Data vector length must be a power of 2");
    }

    int n = a.size();

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

    if (invert) {
        for (cmpx_data_t & x : a)
            x /= n;
    }
}


/**
 * Inherts docs
 */
void fft_recursive(vector<cmpx_data_t> &a, bool invert){

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

    fft_recursive(a0, invert);
    fft_recursive(a1, invert);

    double ang = - 2 * M_PI / n * (invert? -1:1);

    cmpx_data_t w(1), wn(cos(ang), sin(ang));

    for (int i =0; 2*i < n; i++){
        a[i] = a0[i] + w * a1[i];
        a[i + n/2] = a0[i] - w * a1[i];
        if (invert){
            a[i]/= 2;
            a[i + n/2] /= 2;
        }
        w *= wn;
    }
}



/**
 * Inherts docs
 */
void fft2d(vector<vector<cmpx_data_t>> &a, bool invert) {
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
        fft(a[row_index], invert);
    }

    // DO 1D FFT on the resulting array
    for (int col_index = 0; col_index < n; col_index ++){
        // cpy col to a vector
        vector<cmpx_data_t> col;
        for (int i = 0; i<m; i++){
            col.push_back(a[i][col_index]);
        }
        // call fft on col
        fft(col, invert);

        // copy back to matrix
        for (int i = 0; i<m; i++){
            a[i][col_index] = col[i];
        }
    }
}


/**
 * Inherts docs
 */
void fft_data_generator(vector<cmpx_data_t> &input_data, vector<cmpx_data_t> &output_data, vector<float> &frequencies, vector<float> &amplitudes, vector<float> &time_range, bool invert){

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
    fft(output_data, invert);
}



/**
 * Inherts docs
 */
void fft2d_random_data_generator(vector<vector<cmpx_data_t>> &input_data, vector<vector<cmpx_data_t>> &output_data, bool invert){
    int m = input_data.size();
    int n = (input_data[0]).size();

    // Fast Failing Assertions
    if (input_data.size() != output_data.size()){
        throw runtime_error("Input data vector length must match output data vector length");
    }

    if (log2(double(input_data.size())) != int(log2(double(input_data.size())))){
        throw runtime_error("Data vector length must be a power of 2");
    }

    mt19937 gen(1703);
    uniform_real_distribution<> dis(0.0, 10.0);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            double real = dis(gen);
            double imag = dis(gen);
            input_data[i][j] = cmpx_data_t(real, imag);
            output_data[i][j] = cmpx_data_t(real, imag);
        }
    }

    fft2d(output_data, invert);
}


