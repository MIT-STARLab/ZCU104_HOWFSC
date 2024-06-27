/*
 * MIT STAR Lab
 * Modified by Subhi in Jun 13
 */


#include <bits/types/struct_timeval.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/select.h>
#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <time.h>



// XRT includes
#include "xfft_v9_1_bitacc_cmodel.h"
#include "xrt/xrt_bo.h"
#include <experimental/xrt_xclbin.h>
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"


using namespace std;




typedef complex<float> cmpx_data_t;




/**
 * Function to return current time
 */
double get_time_stamp(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_usec + tv.tv_sec * 1e6;
}




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

    if (invert) {
        for (cmpx_data_t & x : a)
            x /= n;
    }
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
 * Checks if two complex numbers are equal within a tolarance
 */
bool approx_equal(const cmpx_data_t &a, const cmpx_data_t &b, double tolerance) {
    return (std::abs(a.real() - b.real()) <= tolerance) &&
           (std::abs(a.imag() - b.imag()) <= tolerance);
}




///////  Testing Data Parameters
const int  FFT_LENGTH = 2048; 
const float sampling_rate = 0.01;
vector<float> frequencies = {1,4,7};
vector<float> amplitudes = {3,1,0.5}; 
vector<float> time_range = {0, sampling_rate * FFT_LENGTH};
bool invert = 0;
double tolerance = 1e-2;


/////// Timing Variables
double hardware_start;
double hardware_end;
double hardware_time;


int main(int argc, char ** argv) {

    cout << "argc = " << argc <<endl;


    ////////////////              Generate Testing Data            /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    cout << "Generating Testing Data..." << endl;
    vector<cmpx_data_t> software_generated_input_data(FFT_LENGTH, 0);      // to include the input data
    vector<cmpx_data_t> software_generated_ouput_data(FFT_LENGTH, 0);      // to include the expected output data

    fft_data_generator(software_generated_input_data, software_generated_ouput_data, frequencies, amplitudes, time_range, invert);

    cout << "Testing Data Generated" << endl;


    ////////////////     Initialize XRT device and load xclbin    //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    int device_index = 0;
    std::string XCLBinaryFileName = argv[1];         // fft_top.xclbin

    std::cout << "Open the device " << device_index << std::endl;
    auto device = xrt::device(device_index);

    std::cout << "Load the xclbin " << XCLBinaryFileName << std::endl;
    auto uuid = device.load_xclbin(XCLBinaryFileName);



    ////////////////        Set Kernel and argument               //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    auto krnl = xrt::kernel(device, uuid, "fft_top", xrt::kernel::cu_access_mode::exclusive);

    size_t input_size_bytes = sizeof(cmpx_data_t) * FFT_LENGTH;

    std::cout << "Allocate Buffer in Global Memory\n";

    auto boIn = xrt::bo(device, input_size_bytes, krnl.group_id(1));
    auto boOut = xrt::bo(device, input_size_bytes, krnl.group_id(2));

    auto boIn_map = boIn.map<cmpx_data_t*>();
    auto boOut_map = boOut.map<cmpx_data_t*>();

    fill(boIn_map, boIn_map + FFT_LENGTH, 0);
    fill(boOut_map, boOut_map + FFT_LENGTH, 0);


    // Fill the data in the buffers
    std::cout << "Filling Buffers with input data\n";
    for (int i = 0; i < FFT_LENGTH; i++) {
        boIn_map[i] = software_generated_input_data[i];
    }
    bool forward_fft = 0;
    bool ovFlo;

    cout << "Synchronize input buffer data to device global memory\n";
    boIn.sync(XCL_BO_SYNC_BO_TO_DEVICE);


    // Execution of the kernel
    cout << "Execution of the kernel\n";
    hardware_start = get_time_stamp();
    auto run = krnl(forward_fft, boIn, boOut, ovFlo);
    run.wait();
    hardware_end = get_time_stamp();
    hardware_time = (hardware_end-hardware_start)/1000;
    std::cout << "Exeution time running 1D FFT with length 2048 on FPGA: " << hardware_time << " msec " << std::endl;



    // Get the output;
    cout << "Get the output data from the device" << std::endl;
    boOut.sync(XCL_BO_SYNC_BO_FROM_DEVICE);


    // Validate Results
    cout << "Validating Output results" << std::endl;
    bool passed = true;

    for (int j = 0; j < FFT_LENGTH; j++) {
        if (!approx_equal(software_generated_ouput_data[j], boOut_map[j], tolerance)) {
            cerr << "Error at index " << j 
                 << ": Expected (" << software_generated_ouput_data[j].real() << ", "  << software_generated_ouput_data[j].imag()
                 << "), Actual (" << boOut_map[j].real() << ", " << boOut_map[j].imag() 
                 << ")" << endl;
            passed = false;
        }
    }


    if (passed) {
        cout << "All tests passed!" << endl;
    } else if (ovFlo){
        cout << " (OVERFLOW!!!)" << endl;
    } else{
        cout << "Tests failed!" << endl;
    }

}