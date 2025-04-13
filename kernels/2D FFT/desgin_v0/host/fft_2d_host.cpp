
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
#include "fft_2d_krnl_include.h"    /* Kernel */
#include "fft.h"                    /* Pure Software Implementation */
// #include "xrt_mem.h"


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




//// MUST DECLARE HERE NOT IN MAIN; OTHERWISE SEG FAULT ISSUES
cmpx_data_t software_generated_input_data[MAT_SIZE];      // to include the input data
cmpx_data_t software_generated_ouput_data[MAT_SIZE];      // to include the expected output data


int main(int argc, char** argv) {

    // Check command line arguments
    if(argc != 2) {
        std::cout << STR_USAGE << argv[0] <<" <binary_container.xclbin> " <<std::endl;
        return EXIT_FAILURE;
    }

    ////////////////              Generate Testing Data            /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    // Compute the size of arrays in bytes
    size_t size_in_bytes = MAT_SIZE * sizeof(cmpx_data_t);
    std::cout << STR_INFO << "Matrix size in words  = " << MAT_SIZE  * 2   << std::endl;
    std::cout << STR_INFO << "Matrix size in bytes  = " << size_in_bytes    << std::endl;


    std::cout << "Generating Testing Data..." << std::endl;
    bool invert = false;
    fft2d_random_data_generator_array(software_generated_input_data, software_generated_ouput_data, MAT_ROWS, MAT_COLS, invert);

    std::cout << "Testing Data Generated Successfully" << std::endl;


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
    auto krnl = xrt::kernel(my_device, xclbin_uuid, "fft_2d:{fft_2d_1}");
    std::cout << STR_PASSED << "auto krnl = xrt::kernel(my_device, xclbin_uuid, \"fft_2d:{fft_2d_1}\")" << std::endl;

    std::cout << STR_INFO << "Allocate Buffers in Global Memory" << std::endl;
                                                                                                  //kernel argument 0: bool direction
    auto bo_input  = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1));      //kernel argument 1: input matrix array
    auto bo_temp   = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2));      //kernel argument 2: temp matrix array
    auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(3));      //kernel argument 3: ouput matrix array


    std::cout << STR_PASSED << "auto bo_input  = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1)) (=" << krnl.group_id(1) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_temp   = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2)) (=" << krnl.group_id(2) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(3)) (=" << krnl.group_id(3) << "))" << std::endl;

    //Map the contents of the buffer object into host memory
    auto bo_input_map  = bo_input.map<cmpx_data_t*>();
    // auto bo_temp_map = bo_temp.map<cmpx_data_t*>(); //DEV ONLY CAN'T MAP
    auto bo_output_map = bo_output.map<cmpx_data_t*>();

    fill(bo_input_map, bo_input_map + MAT_SIZE, cmpx_data_t(0, 0));
    fill(bo_output_map, bo_output_map + MAT_SIZE, cmpx_data_t(0, 0)); 

    std::cout << STR_PASSED << "auto bo_input_map  = bo_input.map<cmpx_data_t*>()" << std::endl;
    std::cout << STR_PASSED << "auto bo_output_map = bo_output.map<cmpx_data_t*>()" << std::endl;

    // Fill the data in the buffers
    std::cout << "Filling Buffers with input data" << std::endl;
    for (int i = 0; i < MAT_SIZE; i++) {
        bo_input_map[i] = software_generated_input_data[i];
    }
    bool forward_fft = 0; //kernel argument 0: direction


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
    auto run = krnl(forward_fft, bo_input, bo_temp, bo_output);
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
    std::cout << "\n************************ TIMING SUMMARY ************************" << std::endl;
    std::cout << "Timing 2D FFT Kernel with input matrix with size (" << MAT_ROWS << ", " << MAT_COLS << ")" << std::endl;
    std::cout << "Time Total = "          << ELAPSED(syncOutputEnd, syncInputTime)    << std::endl;
    std::cout << "Load Input Time  = "    << ELAPSED(syncInputTimeEnd, syncInputTime) << std::endl;
    std::cout << "Excution Time = "       << ELAPSED(runTimeEnd, runTime)             << std::endl;
    std::cout << "Load Output Time = "    << ELAPSED(syncOutputEnd, syncOutputTime)   << std::endl;
    std::cout << "****************************************************************\n" << std::endl;



    ////////////////                 Validation                   //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    cout << "Validating Output results" << std::endl;

    bool passed = true;
    const float tolerance = 1e-2;
    int failed_count = 0;

    float avg_real_relative_err = 0;
    float avg_imag_relative_err = 0;

    for (int j = 0; j < MAT_SIZE; j++) {
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
            avg_real_relative_err +=  (std::abs(( (software_generated_ouput_data[j].real() - bo_output_map[j].real())/ software_generated_ouput_data[j].real())) / MAT_SIZE);
        }

        if (software_generated_ouput_data[j].imag() != 0){
            avg_real_relative_err +=  (std::abs(( (software_generated_ouput_data[j].imag() - bo_output_map[j].imag())/ software_generated_ouput_data[j].imag())) / MAT_SIZE);
        }
    }


    if (passed) {
        cout << "All tests passed!" << endl;
    } else{
        cout << "Tests failed " << "at "<< failed_count << " entries out of " << MAT_SIZE << endl;
    }

    std::cout << "Average Relative Error in the Real Componenet Across All entries = " << avg_real_relative_err << std::endl; 
    std::cout << "Average Relative Error in the Imaginary Componenet Across All entries = " << avg_imag_relative_err << std::endl; 

    return 0;
}