/*
 * MIT STAR Lab
 * M.Subhi Abo Rdan (msubhi_a@mit.edu)
 * Last modified in Dec 29, 2024
 */

#include <sys/syscall.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <complex>
#include <cmath>

// XRT includes
#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

// My Utils includes
#include "fft.h"
#include "angular_spectrum_propagation.h"


using namespace std;
typedef complex<float> cmpx_data_t;

#define MAT_COLS 1024
#define MAT_ROWS 1024
#define MAT_SIZE (MAT_COLS * MAT_ROWS)


#define TIMER(label)  timespec label; syscall(SYS_clock_gettime, CLOCK_MONOTONIC, &label)
#define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0




/////// Error Profiling Helpers ///////
static const char*    STR_PASSED   = "PASSED:   ";
static const char*    STR_FAILED   = "FAILED:   ";
static const char*    STR_INFO     = "INFO:     ";
static const char*    STR_USAGE    = "USAGE:    ";
static const char*    STR_RESULTS  = "RESULTS:  ";

static const char* error_message =
    "Error:\tResults Mismatch:\n"
    "\tIndex i = %d; CPU Results (%f, %f); Device Result (%f, %f)\n"
    "\tRelative error in real component = %f\n"
    "\tRelative error in imag component = %f\n";

static const char* timing_summary = 
    "Round Run Time:\n"
    "-------------------------\n"
    "Matrix Dimension Size =   (%d, %d)\n"
    "Total CPU Time        =   %f\n"
    "Total FPGA Time       =   %f\n"
    "->Input Load Time     =   %f\n"
    "->Excution Time       =   %f\n"
    "->Output Load Time    =   %f\n";

static const char* final_results =
    "\n\n**************************      Timing Summary   ****************************\n"
    "Kernel name                               = Angular Spectrum Propagation\n"
    "Input Data Size                           = (%d, %d)\n"
    "Number of Run Rounds                      = %d\n"
    "FPGA Geometric Mean Run Time (With Sync)  = %f\n"
    "FPGA Geometric Mean Run Time (No   Sync)  = %f\n"
    "CPU  Geometric Mean Run Time              = %f\n"
    "**********************************************************************************\n";


/**
 * @brief Checks if two complex numbers are equal within a relative tolarance
 */
bool approx_equal(const cmpx_data_t &a, const cmpx_data_t &b, double tolerance) {
    bool real_compare = (a.real()<tolerance && b.real()<tolerance) || (abs((a.real() - b.real())/max(abs(a.real()), abs(b.real()))) < tolerance);
    bool imag_compare = (a.imag()<tolerance && b.imag()<tolerance) || (abs((a.imag() - b.imag())/max(abs(a.imag()), abs(b.imag()))) < tolerance);

    return  real_compare && imag_compare;
}

/**
 * Validate results
 */
void validate_results(cmpx_data_t *expected_output, cmpx_data_t *real_output, int data_size, double tolerance, bool save_to_files=false){
    std::cout << STR_INFO << "Start Validation" << std::endl;

    bool passed = true;
    int failed_count = 0;

    float avg_real_relative_err = 0;
    float avg_imag_relative_err = 0;

    for (int j = 0; j < data_size; j++) {

        if (!approx_equal(expected_output[j], real_output[j], tolerance)) {
            printf(error_message, j, expected_output[j].real(), expected_output[j].imag(), real_output[j].real(), real_output[j].imag()
                                , std::abs( (expected_output[j].real()-real_output[j].real()) / (expected_output[j].real()) )
                                , std::abs( (expected_output[j].imag()-real_output[j].imag()) / (expected_output[j].imag()) ));
            passed = false;
            failed_count++;
        }

        if (expected_output[j].real() != 0){
            avg_real_relative_err +=  std::abs( (expected_output[j].real()-real_output[j].real()) / (expected_output[j].real()) ) /data_size;
        }

        if (expected_output[j].imag() != 0){
            avg_imag_relative_err +=  std::abs( (expected_output[j].imag()-real_output[j].imag()) / (expected_output[j].imag()) )/ data_size;
        }
    }

    std::cout << STR_RESULTS << "Average Relative Error in the Real Componenet Across All Indicies = " << avg_real_relative_err << std::endl; 
    std::cout << STR_RESULTS << "Average Relative Error in the Imag Componenet Across All Indicies = " << avg_imag_relative_err << std::endl; 

    if (passed){std::cout << STR_PASSED << "Validation; data at all indicies match";}
    else       {std::cout << STR_FAILED << "Validation; test failed at " << failed_count << "indicies out of " << data_size << std::endl;}

    if (save_to_files){
        ofstream hardware_output_file("hardware_output.txt");
        ofstream software_output_file("software_output.txt");
        hardware_output_file << MAT_ROWS << endl;
        software_output_file << MAT_ROWS << endl;

        for (int i = 0; i < MAT_SIZE; i++) {
            hardware_output_file << hardware_output_data[i].real() << " " << hardware_output_data[i].imag() << endl;
            software_output_file << software_generated_data[i].real() << " " << software_generated_data[i].imag() << endl;
        }
        hardware_output_file.close();
        software_output_file.close();
    }
}



int main(int argc, char** argv) {

    // Check command line arguments
    if(argc < 2) {
        std::cout << STR_USAGE << argv[0] <<" <binary_container.xclbin> <trails>" <<std::endl;
        return EXIT_FAILURE;
    }
    const int rounds = (argc < 3)? 10: std::stoi(argv[2]);

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
    // HLS interface
    // void angular_spectrum(
    //              bool direction,
    //              double scale,
    //              double distance,
    //              double k_2,
    //              double *kxy,
    //              cmpx_data_t *input_mat,
    //              cmpx_data_t *output_mat,
    //              cmpx_data_t *temp_mat_1, 
    //              cmpx_data_t *temp_mat_2
    // );

    // Starting the (HLS) PL kernel
    auto krnl = xrt::kernel(my_device, xclbin_uuid, "angular_spectrum:{angular_spectrum_1}");
    std::cout << STR_PASSED << "auto krnl = xrt::kernel(my_device, xclbin_uuid, \"angular_spectrum:{angular_spectrum_1}\")" << std::endl;

    std::cout << STR_INFO << "Allocate Buffers in Global Memory" << std::endl;

    // Compute the size of arrays in bytes
    size_t size_in_bytes = MAT_SIZE * sizeof(cmpx_data_t);
    std::cout << STR_INFO << "Matrix size in words  = " << MAT_SIZE  * 2   << std::endl;
    std::cout << STR_INFO << "Matrix size in bytes  = " << size_in_bytes   << std::endl;
                                                                                                        // kernel argument 0: bool direction
                                                                                                        // kernel argument 1: double scale
                                                                                                        // kernel argument 2: bool distance
                                                                                                        // kernel argument 3: bool k_2
    auto bo_kxy   = xrt::bo(my_device, MAT_ROWS * sizeof(double), XCL_BO_FLAGS_NONE, krnl.group_id(4)); // kernel argument 4: double *kxy,
    auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(5));             // kernel argument 5: input  matrix array
    auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(6));            // kernel argument 6: output matrix array
    auto bo_temp1 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(7));             // kernel argument 7: temp_mat_1 array
    auto bo_temp2 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(8));             // kernel argument 8: temp_mat_2 array

    std::cout << STR_PASSED << "auto bo_kxy   = xrt::bo(my_device, MAT_ROWS * sizeof(double), XCL_BO_FLAGS_NONE, krnl.group_id(4)) (=" << krnl.group_id(1) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(5)) (=" << krnl.group_id(5) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_ouput = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(6)) (=" << krnl.group_id(6) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_temp1 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(7)) (=" << krnl.group_id(7) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_temp2 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(8)) (=" << krnl.group_id(8) << "))" << std::endl;

    //Map the contents of the buffer object into host memory
    auto bo_kxy_map     = bo_kxy.map<double*>();
    auto bo_input_map   = bo_input.map<cmpx_data_t*>();
    auto bo_output_map  = bo_output.map<cmpx_data_t*>();

    std::cout << STR_PASSED << "auto bo_kxy_map    = bo_kxy.map<double*>()"         << std::endl;
    std::cout << STR_PASSED << "auto bo_input_map  = bo_input.map<cmpx_data_t*>()"  << std::endl;
    std::cout << STR_PASSED << "auto bo_output_map = bo_output.map<cmpx_data_t*>()" << std::endl;



    ////////////////              Generate Testing Data            /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    std::cout << STR_INFO << "Generating Testing Data..." << std::endl;

    bool forward_fft   = 1;
    // input parameters
    double sigma       = 0.2 * MAT_ROWS;
    double intensity   = 1e5;
    double noise_stddev= 0.01;

    // telescope parameters
    double wavelength  = 500e-9;
    double pixel_scale = 10e-3 / (MAT_ROWS / 2);
    double distance    = 1000e-3;

    double delkx = 2.0 * M_PI / (pixel_scale * MAT_ROWS);
    double k = 2.0 * M_PI / wavelength;
    double k_2 = k*k;
    double scale = (1/(double)MAT_SIZE);

    double* kxy = new double[MAT_ROWS];
    cmpx_data_t* software_generated_input_data  = new cmpx_data_t[MAT_SIZE];
    cmpx_data_t* software_generated_output_data = new cmpx_data_t[MAT_SIZE];

    generate_star_gaussian(software_generated_input_data, MAT_ROWS, sigma, intensity, noise_stddev, true);
    for (int i = 0; i < MAT_ROWS; i++){
        kxy[i] = ((-(double)MAT_ROWS / 2.0 + i) + 0.5) * delkx; // Center bins
    }

    std::cout << STR_PASSED << "Testing Data Generated" << std::endl;

    std::cout << STR_INFO << "Filling Argument Buffers with input data" << std::endl;
    for (int i = 0; i < MAT_SIZE; i++) {
        bo_input_map[i] = software_generated_input_data[i];
    }
    for (int i = 0; i < MAT_ROWS; i++) {
        bo_kxy_map[i] = kxy[i];
    }


    ////////////////              Kernel Excution                 //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    std::cout << STR_INFO << "Running on CPU First" << std::endl;
    TIMER(cpuTimeStart);
    angular_spectrum_propagation(software_generated_input_data, MAT_ROWS, (float) wavelength, (float) distance, (float) pixel_scale);
    TIMER(cpuTimeEnd);
    std::cout << STR_PASSED << "CPU Done" << std::endl;

    // Synchronize buffer content with device side
    std::cout << STR_INFO << "synchronize input buffer data to device global memory" << std::endl;
    TIMER(syncInputTime);
    bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_kxy.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    TIMER(syncInputTimeEnd);
    std::cout << STR_PASSED << "bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE)" << std::endl;
    std::cout << STR_PASSED << "bo_kxy.sync(XCL_BO_SYNC_BO_TO_DEVICE)" << std::endl;


    std::cout << STR_INFO <<  "Execution of the kernel" << std::endl;
    TIMER(runTime);
    auto run = krnl(
                forward_fft, 
                scale,
                distance,
                k_2,
                bo_kxy,
                bo_input,
                bo_output,
                bo_temp1, 
                bo_temp2
                );
    std::cout << std::endl << STR_INFO << "Waiting for kernels to end..." << std::endl << std::endl;
    run.wait();
    TIMER(runTimeEnd);
    std::cout << STR_PASSED << "run.wait()" << std::endl;

    // Retrieving the results
    TIMER(syncOutputTime);
    bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    TIMER(syncOutputEnd);
    std::cout << STR_PASSED << "bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE)" << std::endl;


    const double tolerance = 1e-2;
    validate_results(software_generated_output_data, bo_output_map, MAT_SIZE, tolerance);

    // TIMING
    printf(timing_summary, MAT_ROWS, MAT_COLS, 
                           ELAPSED(cpuTimeEnd, cpuTimeStart),
                           ELAPSED(syncOutputEnd, syncInputTime), 
                           ELAPSED(syncInputTimeEnd, syncInputTime), 
                           ELAPSED(runTimeEnd, runTime), 
                           ELAPSED(syncOutputEnd, syncOutputTime));

    return 0;
}
