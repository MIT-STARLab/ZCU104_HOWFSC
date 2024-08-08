
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include "fft_2d_krnl_include.h"

using namespace std;

typedef complex<float> cmpx_data_t;




/**
 * Inherts docs
 */
void fft(vector<cmpx_data_t> &a, bool invert) {

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
void fft2d(vector<vector<cmpx_data_t>> &a, bool invert) {
    int m = a.size();
    int n = (a[0]).size();

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



void fft2d_random_data_generator_array(cmpx_data_t *input_data, cmpx_data_t *output_data, int m, int n, bool invert) {
    mt19937 gen(1703);
    uniform_real_distribution<> dis(100.0, 1000.0);

    vector<vector<cmpx_data_t>> output_data_v(m, vector<cmpx_data_t>(n));

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int index = i * n + j;
            float real = dis(gen);
            float imag = dis(gen);
            input_data[index] = cmpx_data_t(real, imag);
            output_data_v[i][j] = cmpx_data_t(real, imag);
        }
    }

    fft2d(output_data_v, invert);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int index = i * n + j;
            output_data[index] = output_data_v[i][j];
        }
    }
}


/**
 * Checks if two complex numbers are equal within a tolarance
 */
bool approx_equal(const cmpx_data_t &a, const cmpx_data_t &b, double tolerance) {
    bool real_compare = (a.real() == b.real() && b.real()==0) || (std::abs((a.real() - b.real())/a.real()) < tolerance);
    bool imag_compare = (a.imag() == b.imag() && b.imag()==0) || (std::abs((a.imag() - b.imag())/a.imag()) < tolerance);
    return  real_compare && imag_compare;
}

//// MUST DECLARE HERE NOT IN MAIN; OTHERWISE SEG FAULT ISSUES
cmpx_data_t software_generated_input_data[MAT_SIZE];      // to include the input data
cmpx_data_t software_generated_ouput_data[MAT_SIZE];      // to include the expected output data

cmpx_data_t hardware_input_data[MAT_SIZE];
cmpx_data_t hardware_output_data[MAT_SIZE];
cmpx_data_t hardware_temp_data[MAT_SIZE];



int main() {

    cout << "Generate Testing Data using Pure software implementaion" << endl;
    /////////////   Generate Testing Data using Pure software implementaion  //////////////


    fft2d_random_data_generator_array(software_generated_input_data, software_generated_ouput_data, MAT_ROWS, MAT_COLS, false);
    cout << "Testing data generated" << endl;



    cout << "Call FFT function" << endl;
    /////////////////  Call FFT function   //////////////
    bool FORWARD_FFT = true;                                // Use true for forward FFT, false for inverse FFT

    for (int i =0; i<MAT_SIZE; i++){
        hardware_input_data[i] = software_generated_input_data[i];
        hardware_output_data[i] = software_generated_input_data[i];
        hardware_temp_data[i] = software_generated_input_data[i];
    }

    fft_2d(FORWARD_FFT, hardware_input_data, hardware_temp_data, hardware_output_data);



    cout << "Validating Results" << endl;
    /////////////////  Validate Results   //////////////    
    bool passed = true;
    const float tolerance = 1e-2;
    int failed_count = 0;

    float avg_real_relative_err = 0;
    float avg_imag_relative_err = 0;

    for (int j = 0; j < MAT_SIZE; j++) {
        if (!approx_equal(hardware_output_data[j], software_generated_ouput_data[j], tolerance)) {
            passed = false;
            cout << "Error at index " << j 
                 << ": Expected (" << software_generated_ouput_data[j].real() << ", " << software_generated_ouput_data[j].imag() 
                 << "), Actual ("  <<   hardware_output_data[j].real() << ", " <<   hardware_output_data[j].imag() 
                 << ")" << endl;
            cout << "relative real error is " << std::abs(( software_generated_ouput_data[j].real() - hardware_output_data[j].real())/ software_generated_ouput_data[j].real()) << endl;
            cout << "relative imag error is " << std::abs(( software_generated_ouput_data[j].imag() - hardware_output_data[j].imag())/ software_generated_ouput_data[j].imag()) << endl;
            failed_count++;
        }

        if (software_generated_ouput_data[j].real() != 0){
            avg_real_relative_err +=  (std::abs(( software_generated_ouput_data[j].real() - hardware_output_data[j].real()))/software_generated_ouput_data[j].real() / MAT_SIZE);
        }

        if (software_generated_ouput_data[j].imag() != 0){
            avg_real_relative_err +=  (std::abs(( software_generated_ouput_data[j].imag() - hardware_output_data[j].imag())/software_generated_ouput_data[j].imag()) / MAT_SIZE);
        }

    }


    if (passed) {
        cout << "All tests passed!" << endl;
    } else{
        cout <<"Tests failed " << "at "<<failed_count << " indicies out of " << MAT_SIZE << endl << endl;
    }

    
    std::cout << "Average Relative Error in the Real Componenet Across All Indicies = " << avg_real_relative_err << std::endl; 
    std::cout << "Average Relative Error in the Imaginary Componenet Across All Indicies = " << avg_imag_relative_err << std::endl; 



    return 0;
}
