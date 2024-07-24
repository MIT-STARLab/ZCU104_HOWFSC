// Copyright © 2023 Advanced Micro Devices, Inc. All rights reserved.
// SPDX-License-Identifier: MIT

#include <sys/syscall.h>      /* Definition of SYS_* constants */
#include <time.h>

#include <stdlib.h>
#include <iostream>
#include <string>

#include <complex>
#include <cmath>

// XRT includes
#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

// FFT includes
#include "fft_krnl_include.h"    /* Kernel */
#include "fft.h"               /* Pure Software Implementation */


using namespace std;

typedef complex<float> cmpx_data_t;




/**
 * Checks if two complex numbers are equal within a tolarance
 */
bool approx_equal(const cmpx_data_t &a, const cmpx_data_t &b, double tolerance) {
    bool real_compare = (a.real() == b.real() && b.real()==0) || (abs((a.real() - b.real())/a.real()) < tolerance);
    bool imag_compare = (a.imag() == b.imag() && b.imag()==0) || (abs((a.imag() - b.imag())/a.imag()) < tolerance);
    return  real_compare && imag_compare;
}




/////// Timing Helpers ///////
#define TIMER(label)  timespec label; syscall(SYS_clock_gettime, CLOCK_MONOTONIC, &label)
#define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0



/////// Error Profiling Helpers ///////
static const char*    STR_PASSED  = "PASSED:  ";
static const char*    STR_INFO    = "INFO:    ";
static const char*    STR_USAGE   = "USAGE:   ";



///////  Testing Data Parameters  ///////
float sampling_rate = 0.01;
std::vector<float> frequencies = {1,4,7};
std::vector<float> amplitudes = {3,1,0.5}; 
std::vector<float> time_range = {0, sampling_rate * DATA_SIZE};
bool invert = 0;








int main(int argc, char** argv) {

    // Check command line arguments
    if(argc != 2) {
        std::cout << STR_USAGE << argv[0] <<" <binary_container.xclbin> " <<std::endl;
        return EXIT_FAILURE;
    }


    ////////////////              Generate Testing Data            /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    // Compute the size of arrays in bytes
    size_t size_in_bytes = DATA_SIZE * sizeof(cmpx_data_t);
    std::cout << STR_INFO << "DATA size in words  = " << DATA_SIZE  * 2   << std::endl;
    std::cout << STR_INFO << "DATA size in bytes  = " << size_in_bytes    << std::endl;


    std::cout << "Generating Testing Data..." << std::endl;
    std::vector<cmpx_data_t> software_generated_input_data(DATA_SIZE, 0);
    std::vector<cmpx_data_t> software_generated_ouput_data(DATA_SIZE, 0);

    fft_data_generator(software_generated_input_data, software_generated_ouput_data, frequencies, amplitudes, time_range, invert);

    std::cout << "Testing Data Generated" << std::endl;


    ////////////////     Initialize XRT device and load xclbin    //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    // Loading the xclbin
    char* xclbinFilename = argv[1];

    //open the device and return its handle
    unsigned int dev_index = 0;
    auto my_device = xrt::device(dev_index);
    std::cout << STR_PASSED << "auto my_device = xrt::device(" << dev_index << ")" << std::endl;


    //load the xclbin file from disk. This will return the UUID of the xclbin.
    auto xclbin_uuid = my_device.load_xclbin(xclbinFilename);
    std::cout << STR_PASSED << "auto xclbin_uuid = my_device.load_xclbin(" << xclbinFilename << ")" << std::endl;



    ////////////////        Set Kernel and argument               //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    // Starting the (HLS) PL kernel
    auto krnl = xrt::kernel(my_device, xclbin_uuid, "fft_top:{fft_top_1}");
    std::cout << STR_PASSED << "auto krnl = xrt::kernel(my_device, xclbin_uuid, \"fft_top:{fft_top_1}\")" << std::endl;


    std::cout << STR_INFO << "Allocate Buffers in Global Memory" << std::endl;
    //The buffer type is default buffer object (bo) with host buffer and device buffer. The host buffer is allocated and managed by XRT.
    auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1));   //kernel argument 1: input data array
    auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2));  //kernel argument 2: ouput data array
    std::cout << STR_PASSED << "auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1)) (=" << krnl.group_id(1) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2)) (=" << krnl.group_id(2) << "))" << std::endl;

    //Map the contents of the buffer object into host memory
    auto bo_input_map  = bo_input.map<cmpx_data_t*>();
    auto bo_output_map = bo_output.map<cmpx_data_t*>();

    fill(bo_input_map, bo_input_map + DATA_SIZE, cmpx_data_t(0, 0));
    fill(bo_output_map, bo_output_map + DATA_SIZE, cmpx_data_t(0, 0)); 

    std::cout << STR_PASSED << "auto bo_input_map  = bo_input.map<cmpx_data_t*>()" << std::endl;
    std::cout << STR_PASSED << "auto bo_output_map = bo_output.map<cmpx_data_t*>()" << std::endl;
   
    // Fill the data in the buffers
    std::cout << "Filling Buffers with input data" << std::endl;
    for (int i = 0; i < DATA_SIZE; i++) {
        bo_input_map[i] = software_generated_input_data[i];
    }
    bool forward_fft = 0; //kernel argument 0: direction
    bool ovFlo;           //kernel argument 3: ovflow



    ////////////////              Kernel Excution                 //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    // Synchronize buffer content with device side
    std::cout << STR_INFO << "synchronize input buffer data to device global memory" << std::endl;

    TIMER(syncInputTime);
    bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    TIMER(syncInputTimeEnd);
    std::cout << STR_PASSED << "bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE)" << std::endl;


    
    std::cout << STR_INFO <<  "Execution of the kernel" << std::endl;
    TIMER(runTime);
    auto run = krnl(forward_fft, bo_input, bo_output, ovFlo);
    std::cout << std::endl << STR_INFO << "Waiting for kernels to end..." << std::endl << std::endl;
    run.wait();
    TIMER(runTimeEnd);
    std::cout << STR_PASSED << "run.wait()" << std::endl;
    


    // Retrieving the results
    TIMER(syncOutputTime);
    bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    TIMER(syncOutputEnd);
    std::cout << STR_PASSED << "bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE)" << std::endl;

    // TIMING
    std::cout << "Timing 1D FFT Kernel with input data array of size = " << DATA_SIZE << std::endl;
    std::cout << "Time Total = "          << ELAPSED(syncOutputEnd, syncInputTime)    << std::endl;
    std::cout << "Load Input Time  = "    << ELAPSED(syncInputTimeEnd, syncInputTime) << std::endl;
    std::cout << "Excution Time = "       << ELAPSED(runTimeEnd, runTime)             << std::endl;
    std::cout << "Load Output Time = "    << ELAPSED(syncOutputEnd, syncOutputTime)   << std::endl;



    ////////////////                 Validation                   //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    cout << "Validating Output results" << std::endl;

    bool passed = true;
    const float tolerance = 1e-2;
    int failed_count = 0;

    float avg_real_relative_err = 0;
    float avg_imag_relative_err = 0;

    for (int j = 0; j < DATA_SIZE; j++) {
        if (!approx_equal(bo_output_map[j], software_generated_ouput_data[j], tolerance)) {
            passed = false;
            cout << "Error at index " << j 
                 << ": Expected (" << software_generated_ouput_data[j].real() << ", " << software_generated_ouput_data[j].imag() 
                 << "), Actual ("  <<   bo_output_map[j].real() << ", " <<   bo_output_map[j].imag() 
                 << ")" << endl;
            cout << "relative real error is " << std::abs(( software_generated_ouput_data[j].real() - bo_output_map[j].real())/ software_generated_ouput_data[j].real()) << endl;
            cout << "relative imag error is " << std::abs(( software_generated_ouput_data[j].imag() - bo_output_map[j].imag())/ software_generated_ouput_data[j].imag()) << endl;
            failed_count++;
        }

        if (software_generated_ouput_data[j].real() != 0){
            avg_real_relative_err +=  (std::abs(( software_generated_ouput_data[j].real() - bo_output_map[j].real())) / DATA_SIZE);
        }

        if (software_generated_ouput_data[j].imag() != 0){
            avg_real_relative_err +=  (std::abs(( software_generated_ouput_data[j].imag() - bo_output_map[j].imag())) / DATA_SIZE);
        }
    }


    if (passed) {
        cout << "All tests passed!" << endl;
    } else if (ovFlo){
        cout << " (OVERFLOW!!!)" << endl;
    } else{
        cout <<"Tests failed " << "at "<< failed_count << " indicies out of" << DATA_SIZE << endl;
    }


    std::cout << "Average Relative Error in the Real Componenet Across All Indicies = " << avg_real_relative_err << std::endl; 
    std::cout << "Average Relative Error in the Imaginary Componenet Across All Indicies = " << avg_imag_relative_err << std::endl; 


    return 0;
}
