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
#include "xrt/xrt_bo.h"
#include <experimental/xrt_xclbin.h>
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"


using namespace std;


const int  FFT_LENGTH = 2048; 
const bool FORWARD_FFT = true;


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
 * Function to split a string by delimiter and store in vector
 */
vector<string> split(const string &s, char delimiter) {
    vector<string> tokens;
    string token;
    istringstream tokenStream(s);

    while (getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}







/**
 * Checks if two complex numbers are equal within a tolarance
 */
bool approx_equal(const cmpx_data_t &a, const cmpx_data_t &b, double tolerance) {
    return (std::abs(a.real() - b.real()) <= tolerance) &&
           (std::abs(a.imag() - b.imag()) <= tolerance);
}




// timing variables
double hardware_start;
double hardware_end;
double hardware_time;


int main(int argc, char ** argv) {


    std::cout <<"argc = "<<argc <<std::endl;
    std::cout <<argv[0] <<endl;
    std::cout <<argv[1] <<endl;



    ////////////////////    LOAD DATA    /////////////////////
    //////////////////////////////////////////////////////////

    cmpx_data_t input_array[FFT_LENGTH];
    cmpx_data_t expected_output_array[FFT_LENGTH];
    cmpx_data_t actual_output_array[FFT_LENGTH];
    bool ovflo;

    ifstream file("./ftt_tb_data.csv");   

    if (!file.is_open()) {
        cout << "Error opening file!" << endl;
        return 1;
    }

    string line;
    getline(file, line);                                    // Read header line

    for (int i = 0; i < FFT_LENGTH; ++i) {
        if (!getline(file, line)) {
            cerr << "Error reading line " << i + 2 << " from file!" << endl;
            return 1;
        }

        vector<string> tokens = split(line, ',');

        input_array[i] = cmpx_data_t(stof(tokens[0]), stof(tokens[1]));
        expected_output_array[i] = cmpx_data_t(stof(tokens[2]), stof(tokens[3]));
    }

    file.close();
    cout << "Loaded data successfully" << endl;



    ////////////////     Initialize XRT device and load xclbin    //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    // std::string binaryFile = "./fft_top.xclbin";

    int device_index = 0;
    std::string XCLBinaryFileName = argv[1];
    std::cout << "Open the device " << device_index << std::endl;
    auto device = xrt::device(device_index);
    std::cout << "Load the xclbin " << XCLBinaryFileName << std::endl;
    auto uuid = device.load_xclbin(XCLBinaryFileName);

    ////////////////        Set Kernel and argument               //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    auto krnl = xrt::kernel(device, uuid, "fft_top", xrt::kernel::cu_access_mode::exclusive);

    size_t input_size_bytes = sizeof(cmpx_data_t) * FFT_LENGTH;

    std::cout << "Allocate Buffer in Global Memory\n";
    auto boDirection = xrt::bo(device, sizeof(bool), krnl.group_id(0));
    auto boIn = xrt::bo(device, input_size_bytes, krnl.group_id(1));
    auto boOut = xrt::bo(device, input_size_bytes, krnl.group_id(2));
    auto boOvFlo = xrt::bo(device, sizeof(bool), krnl.group_id(3));


    auto boDirection_map = boDirection.map<bool*>();
    auto boIn_map = boIn.map<cmpx_data_t*>();
    auto boOut_map = boOut.map<cmpx_data_t*>();
    auto boOvFlo_map = boOvFlo.map<bool*>();

    // fill(boDirection_map, boDirection_map + 1, 0);
    // fill(boIn_map, boIn_map + FFT_LENGTH, 0);
    // fill(boOut_map, boOut_map + FFT_LENGTH, 0);
    // fill(boOvFlo_map, boOvFlo_map + 1, 0);


    // Fill the data in the buffers
    *boDirection_map = FORWARD_FFT;
    for (int i = 0; i < FFT_LENGTH; i++) {
        boIn_map[i] = input_array[i];
    }


    // Synchronize buffer content with device side
    cout << "Synchronize input buffer data to device global memory\n";
    boDirection.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    boIn.sync(XCL_BO_SYNC_BO_TO_DEVICE);



    // Execution of the kernel
    cout << "Execution of the kernel\n";
    hardware_start = get_time_stamp();
    auto run = krnl(boDirection, boIn, boOut, boOvFlo);
    run.wait();
    hardware_end = get_time_stamp();
    hardware_time = (hardware_end-hardware_start)/1000;
    std::cout << "Exeution time running 1D FFT with length 2048 on FPGA: " << hardware_time << " msec " << std::endl;



    // Get the output;
    cout << "Get the output data from the device" << std::endl;
    boOut.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    boOvFlo.sync(XCL_BO_SYNC_BO_FROM_DEVICE);

    for (int i = 0; i < FFT_LENGTH; i++) {
        actual_output_array[i] = boOut_map[i];
    }
    ovflo = *boOvFlo_map;


    // Validate Results
    bool passed = true;
    double tolerance = 1e-3;

    for (int j = 0; j < FFT_LENGTH; ++j) {
        if (!approx_equal(expected_output_array[j], actual_output_array[j], tolerance)) {
            cerr << "Error at index " << j << ": Expected (" << expected_output_array[j].real() << ", "
                 << expected_output_array[j].imag() << "), Actual (" << actual_output_array[j].real() << ", "
                 << actual_output_array[j].imag() << ")" << endl;
            passed = false;
        }
    }



    if (passed) {
        cout << "All tests passed!" << endl;
        return 0;
    } else {
        cout << "Tests failed!" << endl;
        return 1;
    }

    // ./fft_app fft_top.xclbin
}