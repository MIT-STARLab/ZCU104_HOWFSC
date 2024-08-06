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
#include "fft.h"                 /* Pure Software Implementation */


using namespace std;

typedef complex<float> cmpx_data_t;
#define DATA_SIZE 2048



/////// Timing Helpers ///////
#define TIMER(label)  timespec label; syscall(SYS_clock_gettime, CLOCK_MONOTONIC, &label)
#define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0



/////// Error Profiling Helpers ///////
static const char*    STR_PASSED  = "PASSED:    ";
static const char*    STR_FAILED   = "FAILED: ";
static const char*    STR_INFO    = "INFO:  ";
static const char*    STR_USAGE   = "USAGE: ";
static const char*    STR_RESULTS   = "RESULTS: ";


static const char* error_message =
    "Error:     Results Mismatch:\n"
    "\tIndex i = %d; CPU Results (%f, %f); Device Result (%f, %f)\n"
    "\tRelative error in real component = %f\n"
    "\tRelative error in imag component = %f\n";

static const char* timing_summary = 
    "**************************   Timing Summary **************************\n"
    "Data Size           =   %d\n"
    "Total CPU Time      =   %f\n"
    "Total FPGA Time     =   %f\n"
    "->Input Load Time   =   %f\n"
    "->Excution Time     =   %f\n"
    "->Output Load Time  =   %f\n"
    "***********************************************************************\n";



/**
 * Checks if two complex numbers are equal within a tolarance
 */
bool approx_equal(const cmpx_data_t &a, const cmpx_data_t &b, double tolerance) {

    bool real_compare = (a.real()<tolerance && b.real()<tolerance) || (abs((a.real() - b.real())/max(abs(a.real()), abs(b.real()))) < tolerance);
    bool imag_compare = (a.imag()<tolerance && b.imag()<tolerance) || (abs((a.imag() - b.imag())/max(abs(a.imag()), abs(b.imag()))) < tolerance);

    return  real_compare && imag_compare;
}



/**
 * Validate results
 */
void validate_results(cmpx_data_t *expected_output, cmpx_data_t *real_output, int data_size, double tolerance){

    std::cout << STR_INFO << "Start Validation" << std::endl;

    bool passed = true;
    int failed_count = 0;

    float avg_real_relative_err = 0;
    float avg_imag_relative_err = 0;

    for (int j = 0; j < DATA_SIZE; j++) {

        if (!approx_equal(expected_output[j], real_output[j], tolerance)) {
            printf(error_message, j, expected_output[j].real(), expected_output[j].imag(), real_output[j].real(), real_output[j].imag()
                                , std::abs( (expected_output[j].real()-real_output[j].real()) / (expected_output[j].real()) )
                                , std::abs( (expected_output[j].imag()-real_output[j].imag()) / (expected_output[j].imag()) ));
            passed = false;
            failed_count++;
        }

        if (expected_output[j].real() != 0){
            avg_real_relative_err +=  std::abs( (expected_output[j].real()-real_output[j].real()) / (expected_output[j].real()) ) /DATA_SIZE;
        }

        if (expected_output[j].imag() != 0){
            avg_imag_relative_err +=  std::abs( (expected_output[j].imag()-real_output[j].imag()) / (expected_output[j].imag()) )/ DATA_SIZE;
        }
    }

    std::cout << STR_RESULTS << "Average Relative Error in the Real Componenet Across All Indicies = " << avg_real_relative_err << std::endl; 
    std::cout << STR_RESULTS << "Average Relative Error in the Imag Componenet Across All Indicies = " << avg_imag_relative_err << std::endl; 

    if (passed){std::cout << STR_PASSED << "Validation; data at all indicies match";}
    else       {std::cout << STR_FAILED << "Validation; test failed at " << failed_count << "indicies out of " << DATA_SIZE << std::endl;}

}






///////  First Test Data Parameters  ///////
float sampling_rate = 0.01;
std::vector<float> frequencies = {0.97, 10.44 , 50.44}; // {20 * fund_freq, 250 * fund_freq,  1000 * fund_freq};
std::vector<float> amplitudes = {203,124,34};
std::vector<float> time_range = {0, sampling_rate * DATA_SIZE};
bool invert = 0;
bool FORWARD_FFT = true;
bool scale = false;  //hardware doesn't scale
const double tolerance = 1e-2;




int main(int argc, char** argv) {


    // Check command line arguments
    if(argc < 2) {
        std::cout << STR_USAGE << argv[0] <<" <binary_container.xclbin> <trails>" <<std::endl;
        return EXIT_FAILURE;
    }


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


    size_t size_in_bytes = DATA_SIZE * sizeof(cmpx_data_t);
    std::cout << STR_INFO << "DATA size in words  = " << size_in_bytes/4  << std::endl;
    std::cout << STR_INFO << "DATA size in bytes  = " << size_in_bytes    << std::endl;


    //The buffer type is default buffer object (bo) with host buffer and device buffer. The host buffer is allocated and managed by XRT.
    auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1));   //kernel argument 1: input data array
    auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2));  //kernel argument 2: ouput data array
    std::cout << STR_PASSED << "auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1)) (=" << krnl.group_id(1) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2)) (=" << krnl.group_id(2) << "))" << std::endl;


    //Map the contents of the buffer object into host memory
    auto bo_input_map  = bo_input.map<cmpx_data_t*>();
    auto bo_output_map = bo_output.map<cmpx_data_t*>();

    std::cout << STR_PASSED << "auto bo_input_map  = bo_input.map<cmpx_data_t*>()" << std::endl;
    std::cout << STR_PASSED << "auto bo_output_map = bo_output.map<cmpx_data_t*>()" << std::endl;





    ///////////////////        Generate First Testing Data         /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    std::cout << STR_INFO << "Generating Testing Data..." << std::endl;

    cmpx_data_t* software_generated_input_data = new cmpx_data_t[DATA_SIZE];
    cmpx_data_t* software_generated_output_data = new cmpx_data_t[DATA_SIZE];
    fft_data_generator(software_generated_input_data, software_generated_output_data, DATA_SIZE, frequencies, amplitudes, time_range, invert, scale);

    std::cout << STR_PASSED << "Testing Data Generated" << std::endl;

    std::cout << STR_INFO << "Filling Argument Buffers with input data" << std::endl;

    fill(bo_input_map, bo_input_map + DATA_SIZE, cmpx_data_t(0, 0));
    fill(bo_output_map, bo_output_map + DATA_SIZE, cmpx_data_t(0, 0)); 

    for (int i = 0; i < DATA_SIZE; i++) {
        bo_input_map[i] = software_generated_input_data[i];
    }
    bool forward_fft = FORWARD_FFT;     //kernel argument 0: direction





    ///////////////////////////   First Kernel Excution        /////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    std::cout << STR_INFO << "Running on CPU First" << std::endl;
    TIMER(cpuTimeStart);
    fft(software_generated_input_data, DATA_SIZE, invert, scale);
    TIMER(cpuTimeEnd);
    std::cout << STR_PASSED << "CPU Done" << std::endl;

    // Synchronize buffer content with device side
    std::cout << STR_INFO << "synchronize input buffer data to device global memory" << std::endl;

    TIMER(syncInputTime);
    bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    TIMER(syncInputTimeEnd);
    std::cout << STR_PASSED << "bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE)" << std::endl;


    std::cout << STR_INFO <<  "Execution of the kernel" << std::endl;
    TIMER(runTime);
    auto run = krnl(forward_fft, bo_input, bo_output);
    std::cout << std::endl << STR_INFO << "Waiting for kernels to end..." << std::endl << std::endl;
    run.wait();
    TIMER(runTimeEnd);
    std::cout << STR_PASSED << "run.wait()" << std::endl;


    // Retrieving the results
    TIMER(syncOutputTime);
    bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    TIMER(syncOutputEnd);
    std::cout << STR_PASSED << "bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE)" << std::endl;

    
    validate_results(software_generated_output_data, bo_output_map, DATA_SIZE, tolerance);

    // TIMING
    printf(timing_summary, DATA_SIZE, 
                           ELAPSED(cpuTimeEnd, cpuTimeStart),
                           ELAPSED(syncOutputEnd, syncInputTime), 
                           ELAPSED(syncInputTimeEnd, syncInputTime), 
                           ELAPSED(runTimeEnd, runTime), 
                           ELAPSED(syncOutputEnd, syncOutputTime));


    ///////////////     Multiple Kernel Runs on Random Data for timing      ////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    const int rounds = (argc < 3)? 10: std::stoi(argv[2]);
    std::cout << STR_INFO << "RUN KERNEL NUMBER OF TIMES =  " << rounds << std::endl;

    double device_time_average = 0;
    double cpu_time_average = 0;

    for (int r=0; r<rounds; r++){
        std::cout << STR_INFO << "STATING ROUND : " << r << std::endl;

        fft_random_data_generator(software_generated_input_data, software_generated_output_data, DATA_SIZE, invert, scale = false);

        for (int i = 0; i < DATA_SIZE; i++) {
            bo_input_map[i] = software_generated_input_data[i];
        }

        TIMER(cpuTimeStart);
        fft(software_generated_input_data, DATA_SIZE, invert, scale);
        TIMER(cpuTimeEnd);

        cpu_time_average += ELAPSED(cpuTimeEnd, cpuTimeStart)/rounds;


        TIMER(syncInputTime);
        bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        TIMER(syncInputTimeEnd);

        TIMER(runTime);
        run.start();
        run.wait();
        TIMER(runTimeEnd);

        TIMER(syncOutputTime);
        bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
        TIMER(syncOutputEnd);

        validate_results(software_generated_output_data, bo_output_map, DATA_SIZE, tolerance);

        printf(timing_summary,  DATA_SIZE, 
                                ELAPSED(cpuTimeEnd, cpuTimeStart),
                                ELAPSED(syncOutputEnd, syncInputTime), 
                                ELAPSED(syncInputTimeEnd, syncInputTime), 
                                ELAPSED(runTimeEnd, runTime), 
                                ELAPSED(syncOutputEnd, syncOutputTime));
        
        device_time_average += ELAPSED(syncOutputEnd, syncInputTime)/rounds;
    }

    std::cout << STR_PASSED << "RUN KERNEL NUMBER OF TIMES =  " << rounds << std::endl;
    std::cout << STR_RESULTS << "FPGA Time average  = " << device_time_average << std::endl;
    std::cout << STR_RESULTS << "CPU  Time average  = " << cpu_time_average << std::endl;


    delete[] software_generated_input_data;
    delete[] software_generated_output_data;

    return 0;
}