/*
 * MIT STAR Lab
 * Nicholas Belsten (nbelsten@mit.edu)
 * Last modified on February 2025
 * Host for QR Decomposition Kernel on ZCU104 MPSoC
 */

#include <cstring>
#include <sys/syscall.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

// XRT includes
#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

// Utils includes
#include "qrdcmp.h"





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
    "Round %d Run Time:\n"
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

    "FPGA Arithmetic Mean Run Time (With Sync)  = %f\n"
    "FPGA Arithmetic Mean Run Time (No   Sync)  = %f\n"
    "CPU  Arithmetic Mean Run Time              = %f\n"

    "**********************************************************************************\n";


/**
 * Checks if two float values are approximately equal within a tolerance.
 */
bool approx_equal(float a, float b, float tolerance) {
    if (std::abs(a) < tolerance && std::abs(b) < tolerance) {
        // If both values are near zero, consider them approximately equal
        return true;
    }

    // Compare the relative difference
    return std::abs((a - b) / std::max(std::abs(a), std::abs(b))) < tolerance;
}

/**
 * Validate results for real float arrays.
 */
void validate_results(float *expected_output, float *real_output, int data_size, float tolerance, bool save_to_files, bool print_errors) {
    std::cout << STR_INFO << "Start Validation" << std::endl;

    bool passed = true;
    int failed_count = 0;
    float avg_relative_error = 0;

    for (int j = 0; j < data_size; j++) {
        if (!approx_equal(expected_output[j], real_output[j], tolerance)) {
            if (print_errors) {
                printf(error_message, j, expected_output[j], real_output[j], 
                       std::abs((expected_output[j] - real_output[j]) / expected_output[j]));
            }
            passed = false;
            failed_count++;
        }

        if (expected_output[j] > tolerance) {
            avg_relative_error += std::abs((expected_output[j] - real_output[j]) / expected_output[j]) / data_size;
        }
    }

    std::cout << STR_RESULTS << "Average Relative Error Across All Indices = " << avg_relative_error << std::endl;

    if (passed) {
        std::cout << STR_PASSED << "Validation; data at all indices match" << std::endl;
    } else {
        std::cout << STR_FAILED << "Validation; test failed at " << failed_count << " indices out of " << data_size << std::endl;
    }

    if (save_to_files) {
        std::ofstream hardware_output_file("hardware_output.txt");
        std::ofstream software_output_file("software_output.txt");

        hardware_output_file << N << std::endl;
        software_output_file << N << std::endl;

        for (int i = 0; i < data_size; i++) {
            hardware_output_file << real_output[i] << std::endl;
            software_output_file << expected_output[i] << std::endl;
        }

        hardware_output_file.close();
        software_output_file.close();
    }
}


int main(int argc, char** argv) {

    // Check command line arguments
    if(argc < 2) {
        std::cout << STR_USAGE << argv[0] <<" <binary_container.xclbin> <trials> <print_errors>" <<std::endl;
        return EXIT_FAILURE;
    }
    int rounds = (argc < 3)? 1: std::stoi(argv[2]);
    bool print_errors = (argc < 4)? 0: std::stoi(argv[3]);

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
    //  void krnl_qr_dcmp(data_t *Q, data_t *R);
    // 

    // Starting the (HLS) PL kernel
    auto krnl = xrt::kernel(my_device, xclbin_uuid, "krnl_qr_dcmp:{krnl_qr_dcmp_1}");
    std::cout << STR_PASSED << "auto krnl = xrt::kernel(my_device, xclbin_uuid, \"krnl_qr_dcmp:{krnl_qr_dcmp_1}\")" << std::endl;

    std::cout << STR_INFO << "Allocate Buffers in Global Memory" << std::endl;

    // Compute the size of arrays in bytes
    size_t size_in_bytes = MAT_SIZE * sizeof(data_t);
    std::cout << STR_INFO << "Matrix size in words  = " << MAT_SIZE  * 2   << std::endl;
    std::cout << STR_INFO << "Matrix size in bytes  = " << size_in_bytes   << std::endl;

    auto bo_QT = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(0));             // kernel argument 0: input  matrix array
    auto bo_R = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1));            // kernel argument 1: output matrix array

    std::cout << STR_PASSED << "auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(0)) (=" << krnl.group_id(0) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_ouput = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1)) (=" << krnl.group_id(1) << "))" << std::endl;

    //Map the contents of the buffer object into host memory
    auto bo_QT_map  = bo_QT.map<data_t*>();
    auto bo_R_map  = bo_R.map<data_t*>();

    std::cout << STR_PASSED << "auto bo_QT_map  = bo_QT.map<data_t*>()"  << std::endl;
    std::cout << STR_PASSED << "auto bo_R_map  = bo_R.map<data_t*>()" << std::endl;



    ////////////////              Generate Testing Data            /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    std::cout << STR_INFO << "Generating Testing Data..." << std::endl;


    data_t* software_A = new data_t[MAT_SIZE];
    data_t* software_QT = new data_t[MAT_SIZE];
    data_t* software_R = new data_t[MAT_SIZE];

    data_t* hardware_A = new data_t[MAT_SIZE];
    data_t* hardware_QT = new data_t[MAT_SIZE];
    data_t* hardware_R = new data_t[MAT_SIZE];

    cout << "Generate Testing Data using Pure software implementaion" << endl;
    /////////////   Generate Testing Data using Pure software implementaion  //////////////


    random_data_generator_array(software_A);
    // for(int i=0; i<MAT_SIZE; i++){
    //     software_A[i]=i+1;
    // }
    eye(software_QT);

    for(int i=0; i<MAT_SIZE; i++){
        software_R[i]=software_A[i];

        hardware_A[i]=software_A[i];
        hardware_R[i]=software_R[i];
        hardware_QT[i]=software_QT[i];
   }

    /////////////////  For software QR DCMP function   //////////////
    Matrix A, QT, R;
    A.data = software_A;
    A.cols = N;
    A.rows = N;

    QT.data = software_QT;
    QT.cols = N;
    QT.rows = N;

    R.data=software_R;
    R.cols=N;
    R.rows=N;



    std::cout << STR_PASSED << "Testing Data Generated" << std::endl;




    std::cout << STR_INFO << "Filling Argument Buffers with input data" << std::endl;
    std::memcpy(bo_QT_map, hardware_QT, MAT_SIZE * sizeof(data_t));
    std::memcpy(bo_R_map, hardware_R, MAT_SIZE * sizeof(data_t));



    ////////////////              Kernel Excution                 //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    std::cout << STR_INFO << "Running on CPU First" << std::endl;
    TIMER(cpuTimeStart);
    qr_givens(&A, &QT, &R);
    TIMER(cpuTimeEnd);
    std::cout << STR_PASSED << "CPU Done" << std::endl;

    // Synchronize buffer content with device side
    std::cout << STR_INFO << "synchronize input buffer data to device global memory" << std::endl;
    TIMER(syncInputTime);
    bo_QT.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_R.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    TIMER(syncInputTimeEnd);
    std::cout << STR_PASSED << "bo_QT/R.sync(XCL_BO_SYNC_BO_TO_DEVICE)" << std::endl;

    std::cout << STR_INFO <<  "Execution of the kernel" << std::endl;
    TIMER(runTime);
    auto run = krnl(bo_QT, bo_R);
    run.wait();
    TIMER(runTimeEnd);
    std::cout << STR_PASSED << "run.wait()" << std::endl;

    // Retrieving the results
    TIMER(syncOutputTime);
    bo_QT.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    bo_R.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    TIMER(syncOutputEnd);

    std::cout << STR_PASSED << "bo_QT/R.sync(XCL_BO_SYNC_BO_FROM_DEVICE)" << std::endl;

    // TIMING
    printf(timing_summary, 0, N, N, 
                           ELAPSED(cpuTimeEnd, cpuTimeStart),
                           ELAPSED(syncOutputEnd, syncInputTime), 
                           ELAPSED(syncInputTimeEnd, syncInputTime), 
                           ELAPSED(runTimeEnd, runTime), 
                           ELAPSED(syncOutputEnd, syncOutputTime));

    double tolerance = 1e-5*N;

    std::cout << STR_INFO << "Validating QT" << std::endl;
    validate_results(software_QT, bo_QT_map, MAT_SIZE, tolerance, true, print_errors);
    std::cout << STR_INFO << "Validating R" << std::endl;
    validate_results(software_R, bo_R_map, MAT_SIZE, tolerance, true, print_errors);



    ///////////////     Multiple Kernel Runs on Data to average timing      ////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    if (rounds == 0) {return 0;};

    // std::cout << STR_INFO << "RUN KERNEL NUMBER OF TIMES =  " << rounds << std::endl;

    // double device_exc_time_arithmetic_mean = 0;
    // double device_tot_time_arithmetic_mean = 0;
    // double cpu_time_arithmetic_mean = 0;

    // double device_exc_time_geometric_mean = 1;
    // double device_tot_time_geometric_mean = 1;
    // double cpu_time_geometric_mean = 1;


    // for (int r=0; r<rounds; r++){
    //     std::cout << STR_INFO << "STARTING ROUND : " << r << std::endl;

    //     // change testing parameters each round: maybe?
    //     sigma       = input_sigma;
    //     intensity   = 1e5;
    //     noise_stddev= 0.01;

    //     wavelength  = 500e-9;
    //     pixel_scale = 10e-3 / (MAT_ROWS / 2);
    //     distance    = 1000e-3;
    //     delkx = 2.0 * M_PI / (pixel_scale * MAT_ROWS);
    //     k = 2.0 * M_PI / wavelength;
    //     k_2 = k*k;

    //     generate_star_gaussian(software_generated_input_data, MAT_ROWS, sigma, intensity, noise_stddev);

    //     std::memcpy(bo_input_map, software_generated_input_data, MAT_SIZE * sizeof(cmpx_data_t));

    //     TIMER(cpuTimeStart);
    //     angular_spectrum_propagation(software_generated_input_data, MAT_ROWS, wavelength, distance, pixel_scale);
    //     TIMER(cpuTimeEnd);

    //     cpu_time_arithmetic_mean += ELAPSED(cpuTimeEnd, cpuTimeStart);
    //     cpu_time_geometric_mean *= ELAPSED(cpuTimeEnd, cpuTimeStart);


    //     // propagate any argument changes to kerenl
    //     run.set_arg(1, distance);
    //     run.set_arg(2, k_2);
    //     run.set_arg(3, delkx);

    //     TIMER(syncInputTime);
    //     bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    //     TIMER(syncInputTimeEnd);

    //     TIMER(runTime);
    //     run.start();
    //     run.wait();
    //     TIMER(runTimeEnd);

    //     TIMER(syncOutputTime);
    //     bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    //     TIMER(syncOutputTimeEnd);

    //     device_exc_time_arithmetic_mean += ELAPSED(runTimeEnd, runTime);
    //     device_exc_time_geometric_mean  *= ELAPSED(runTimeEnd, runTime);

    //     device_tot_time_arithmetic_mean += (ELAPSED(syncOutputTimeEnd, syncOutputTime) + ELAPSED(syncInputTimeEnd, syncInputTime) + ELAPSED(runTimeEnd, runTime));
    //     device_tot_time_geometric_mean *= (ELAPSED(syncOutputTimeEnd, syncOutputTime) + ELAPSED(syncInputTimeEnd, syncInputTime) + ELAPSED(runTimeEnd, runTime));

    //     printf(timing_summary,  r+1, MAT_ROWS, MAT_COLS, 
    //                             ELAPSED(cpuTimeEnd, cpuTimeStart),
    //                             ELAPSED(syncOutputTimeEnd, syncOutputTime), 
    //                             ELAPSED(syncInputTimeEnd, syncInputTime), 
    //                             ELAPSED(runTimeEnd, runTime),
    //                             ELAPSED(runTimeEnd, runTime) + ELAPSED(syncOutputTimeEnd, syncOutputTime) + ELAPSED(syncInputTimeEnd, syncInputTime));

    //     validate_results(software_generated_input_data, bo_output_map, MAT_SIZE, tolerance, false, print_errors);
    // }

    // device_exc_time_arithmetic_mean /= rounds;
    // device_tot_time_arithmetic_mean /= rounds;
    // cpu_time_arithmetic_mean /= rounds;

    // device_exc_time_geometric_mean = std::pow(device_exc_time_geometric_mean, 1.0 / rounds);
    // device_tot_time_geometric_mean = std::pow(device_tot_time_geometric_mean, 1.0 / rounds);;
    // cpu_time_geometric_mean = std::pow(cpu_time_geometric_mean, 1.0 / rounds);;

    // printf(final_results, MAT_ROWS, MAT_COLS,
    //                       rounds,

    //                       device_tot_time_geometric_mean,
    //                       device_exc_time_geometric_mean,
    //                       cpu_time_geometric_mean,

    //                       device_tot_time_arithmetic_mean,
    //                       device_exc_time_arithmetic_mean,
    //                       cpu_time_arithmetic_mean);

    return 0;
}
