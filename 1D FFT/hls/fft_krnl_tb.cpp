

/*
 * MIT STAR Lab
 * Written by Subhi in Jun 24
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include "fft_krnl_include.h"


using namespace std;

typedef complex<float> cmpx_data_t;

////// Helpers ///////

/**
 * Implementation of the 1D Cooley-Tukey FFT, the input should have a length of power of 2. Non Recursive Implementation.
 * Implementation taken from https://cp-algorithms.com/algebra/fft.html (accessed on Jun 25).
 * @param a  input vector
 * @param invert 1 if inverse fft
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

    // The hardware FFT IP doesn't scale for inverse also
    // if (invert) {
    //     for (cmpx_data_t & x : a)
    //         x /= n;
    // }
}



/**
 * Produces FFT Testing Data and expected output
 * @param input_data   a vector that will later contatins the input data, must be a power of 2 size
 * @param output_data  a vector that will later contatins the results of applying fft on the input data, must match the size of input_data
 * @param frequencies  a vector containing the frequencies to be included
 * @param amplitudes   a vector containing the corresponding amplitudes, must match the size of frequencies
 * @param time_range   it's first element contatins the timing of first measurment. the second the last measurment.
 * @param invert       whether you want the output inversefft or fft
 */
void fft_data_generator(vector<cmpx_data_t> &input_data, vector<cmpx_data_t> &output_data,const vector<float> &frequencies,const vector<float> &amplitudes, const vector<float> &time_range, bool invert){

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
 * Checks if two complex numbers are equal within a tolarance
 */
bool approx_equal(const cmpx_data_t &a, const cmpx_data_t &b, double tolerance) {
    bool real_compare = (a.real() == b.real() && b.real()==0) || (std::abs((a.real() - b.real())/a.real()) < tolerance);
    bool imag_compare = (a.imag() == b.imag() && b.imag()==0) || (std::abs((a.imag() - b.imag())/a.imag()) < tolerance);
    return  real_compare && imag_compare;
}



/**
 * Writes complex data vector to files
 */
void write_to_file(const vector<cmpx_data_t>& data, const string& filename) {
    ofstream outfile(filename);
    for (const auto& val : data) {
        outfile << val.real() << " " << val.imag() << endl;
    }
    outfile.close();
}




///////  Testing Data Parameters  ///////
const float sampling_rate = 0.01;
const vector<float> frequencies = {0.97, 10.44 , 50.44}; // {20 * fund_freq, 250 * fund_freq,  1000 * fund_freq};
const vector<float> amplitudes = {203,124,34};
const vector<float> time_range = {0, sampling_rate * DATA_SIZE};
bool invert = 1;



int main() {

    cout << "Generate Testing Data using Pure software implementaion" << endl;
    /////////////   Generate Testing Data using Pure software implementaion  //////////////
    vector<cmpx_data_t> software_generated_input_data(DATA_SIZE, 0);      // to include the input data
    vector<cmpx_data_t> software_generated_ouput_data(DATA_SIZE, 0);      // to include the expected output data
    fft_data_generator(software_generated_input_data, software_generated_ouput_data, frequencies,amplitudes, time_range,invert);

    cout << "Saving Testing Data to files software_input.txt, software_ouput.txt" << endl;
    write_to_file(software_generated_input_data, "software_input.txt");
    write_to_file(software_generated_ouput_data, "software_output.txt");
    cout << "Testing data generated" << endl;


    cout << "Call FFT function" << endl;
    /////////////////  Call FFT function   //////////////
    cmpx_data_t hardware_input_data[DATA_SIZE];
    cmpx_data_t hardware_output_data[DATA_SIZE];
    bool ovflo = false;
    bool FORWARD_FFT = false;                                // Use true for forward FFT, false for inverse FFT

    for (int i =0; i<DATA_SIZE; i++){
        hardware_input_data[i] = software_generated_input_data[i];
    }

    fft_top(FORWARD_FFT, hardware_input_data, hardware_output_data, ovflo);

    ofstream hardware_output_file("hardware_output.txt");
    for (int i = 0; i < DATA_SIZE; i++) {
        hardware_output_file << hardware_output_data[i].real() << " " << hardware_output_data[i].imag() << endl;
    }
    hardware_output_file.close();


    cout << "Validating Results" << endl;
    /////////////////  Validate Results   //////////////    
    bool passed = true;
    const float tolerance = 1e-2;
    int failed_count = 0;

    for (int j = 0; j < DATA_SIZE; j++) {
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
    }

    if (passed) {
        cout << "All tests passed!" << endl;
        return 0;
    } else if (ovflo){
        cout << " (OVERFLOW!!!)" << endl;
        return 1;
    } else{
        cout <<"Tests failed " << "at "<<failed_count << " indicies out of" << DATA_SIZE << endl << endl;
        return 0;
    }

}