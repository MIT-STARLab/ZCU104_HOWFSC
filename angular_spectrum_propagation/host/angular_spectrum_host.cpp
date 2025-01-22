/*
 * MIT STAR Lab
 * M.Subhi Abo Rdan (msubhi_a@mit.edu)
 * Last modified on January 16, 2024
 * Pure software implementations of Angular Spectrum Propagation Method to help test hardware kernels
 */

#include <cstring>
#include <sys/syscall.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <cmath>

// XRT includes
#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

// Utils includes
#include "fft.h"
#include "angular_spectrum_propagation.h"


// using namespace std;
typedef float data_t;
typedef std::complex<float> cmpx_data_t;

#define MAT_ROWS 256
#define MAT_COLS 256
#define MAT_SIZE (MAT_ROWS * MAT_COLS)

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

    "FPGA Arithmetic Mean Run Time (With Sync)  = %f\n"
    "FPGA Arithmetic Mean Run Time (No   Sync)  = %f\n"
    "CPU  Arithmetic Mean Run Time              = %f\n"

    "**********************************************************************************\n";


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
void validate_results(cmpx_data_t *expected_output, cmpx_data_t *real_output, int data_size, double tolerance, bool save_to_files, bool print_errors){
    std::cout << STR_INFO << "Start Validation" << std::endl;

    bool passed = true;
    int failed_count = 0;

    float avg_real_relative_err = 0;
    float avg_imag_relative_err = 0;

    for (int j = 0; j < data_size; j++) {

        if (!approx_equal(expected_output[j], real_output[j], tolerance)) {
            if (print_errors){
                printf(error_message, j, expected_output[j].real(), expected_output[j].imag(), real_output[j].real(), real_output[j].imag()
                                , std::abs( (expected_output[j].real()-real_output[j].real()) / (expected_output[j].real()) )
                                , std::abs( (expected_output[j].imag()-real_output[j].imag()) / (expected_output[j].imag()) ));
            }
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
            hardware_output_file << real_output[i].real() << " " << real_output[i].imag() << endl;
            software_output_file << expected_output[i].real() << " " << expected_output[i].imag() << endl;
        }
        hardware_output_file.close();
        software_output_file.close();
    }
}



int main(int argc, char** argv) {

    // Check command line arguments
    if(argc < 2) {
        std::cout << STR_USAGE << argv[0] <<" <binary_container.xclbin> <trails> <print_errors>" <<std::endl;
        return EXIT_FAILURE;
    }
    int rounds = (argc < 3)? 10: std::stoi(argv[2]);
    bool print_errors = (argc < 4 || argc < 3)? 1: std::stoi(argv[3]);

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
    //              data_t scale,
    //              data_t distance,
    //              data_t k_2,
    //              data_t *kxy,
    //              cmpx_data_t *input_mat,
    //              cmpx_data_t *output_mat,
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
                                                                                                        // kernel argument 1: data_t scale
                                                                                                        // kernel argument 2: data_t distance
                                                                                                        // kernel argument 3: data_t k_2
    auto bo_kxy   = xrt::bo(my_device, MAT_ROWS * sizeof(data_t), XCL_BO_FLAGS_NONE, krnl.group_id(4)); // kernel argument 4: data_t *kxy,
    auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(5));             // kernel argument 5: input  matrix array
    auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(6));            // kernel argument 6: output matrix array

    std::cout << STR_PASSED << "auto bo_kxy   = xrt::bo(my_device, MAT_ROWS * sizeof(data_t), XCL_BO_FLAGS_NONE, krnl.group_id(4)) (=" << krnl.group_id(1) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(5)) (=" << krnl.group_id(5) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_ouput = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(6)) (=" << krnl.group_id(6) << "))" << std::endl;

    //Map the contents of the buffer object into host memory
    auto bo_kxy_map     = bo_kxy.map<data_t*>();
    auto bo_input_map   = bo_input.map<cmpx_data_t*>();
    auto bo_output_map  = bo_output.map<cmpx_data_t*>();

    std::cout << STR_PASSED << "auto bo_kxy_map    = bo_kxy.map<data_t*>()"         << std::endl;
    std::cout << STR_PASSED << "auto bo_input_map  = bo_input.map<cmpx_data_t*>()"  << std::endl;
    std::cout << STR_PASSED << "auto bo_output_map = bo_output.map<cmpx_data_t*>()" << std::endl;



    ////////////////              Generate Testing Data            /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    std::cout << STR_INFO << "Generating Testing Data..." << std::endl;

    bool inverse_fft   = false;
    data_t sigma       = 1;
    data_t intensity   = 1e5;
    data_t noise_stddev= 0.01;
    data_t wavelength  = 500e-9;
    data_t pixel_scale = 10e-3 / (MAT_ROWS / 2);
    data_t distance    = 1000e-3;
    data_t delkx = 2.0 * M_PI / (pixel_scale * MAT_ROWS);
    data_t k = 2.0 * M_PI / wavelength;
    data_t k_2 = k*k;
    data_t scale = (1/(data_t)MAT_SIZE);

    data_t* kxy = new data_t[MAT_ROWS];
    cmpx_data_t* software_generated_input_data  = new cmpx_data_t[MAT_SIZE];

    generate_star_gaussian(software_generated_input_data, MAT_ROWS, sigma, intensity, noise_stddev);
    for (int i = 0; i < MAT_ROWS; i++){
        kxy[i] = ((-(data_t)MAT_ROWS / 2.0 + i) + 0.5) * delkx; // Center bins
    }

    std::cout << STR_PASSED << "Testing Data Generated" << std::endl;

    std::cout << STR_INFO << "Filling Argument Buffers with input data" << std::endl;
    std::memcpy(bo_input_map, software_generated_input_data, MAT_SIZE * sizeof(cmpx_data_t));
    std::memcpy(bo_kxy_map, kxy, MAT_ROWS * sizeof(data_t));


    ////////////////              Kernel Excution                 //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    std::cout << STR_INFO << "Running on CPU First" << std::endl;
    TIMER(cpuTimeStart);
    angular_spectrum_propagation(software_generated_input_data, MAT_ROWS, wavelength, distance, pixel_scale);
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
                inverse_fft, 
                scale,
                distance,
                k_2,
                bo_kxy,
                bo_input,
                bo_output
                );
    run.wait();
    TIMER(runTimeEnd);
    std::cout << STR_PASSED << "run.wait()" << std::endl;

    // Retrieving the results
    TIMER(syncOutputTime);
    bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    TIMER(syncOutputEnd);
    std::cout << STR_PASSED << "bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE)" << std::endl;

    // TIMING
    printf(timing_summary, MAT_ROWS, MAT_COLS, 
                           ELAPSED(cpuTimeEnd, cpuTimeStart),
                           ELAPSED(syncOutputEnd, syncInputTime), 
                           ELAPSED(syncInputTimeEnd, syncInputTime), 
                           ELAPSED(runTimeEnd, runTime), 
                           ELAPSED(syncOutputEnd, syncOutputTime));

    double tolerance = 1e-2;
    validate_results(software_generated_input_data, bo_output_map, MAT_SIZE, tolerance, true, print_errors);




    ///////////////     Multiple Kernel Runs on Data to average timing      ////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    if (rounds == 0) {return 0;};

    std::cout << STR_INFO << "RUN KERNEL NUMBER OF TIMES =  " << rounds << std::endl;

    double device_exc_time_arithmetic_mean = 0;
    double device_tot_time_arithmetic_mean = 0;
    double cpu_time_arithmetic_mean = 0;

    double device_exc_time_geometric_mean = 1;
    double device_tot_time_geometric_mean = 1;
    double cpu_time_geometric_mean = 1;


    for (int r=0; r<rounds; r++){
        std::cout << STR_INFO << "STARTING ROUND : " << r << std::endl;

        // change testing parameters each round
        intensity   = 1e5;
        noise_stddev= 0.01;
        wavelength  = 500e-9;
        pixel_scale = 10e-3 / (MAT_ROWS / 2);
        distance    = 1000e-3;
        delkx = 2.0 * M_PI / (pixel_scale * MAT_ROWS);
        k = 2.0 * M_PI / wavelength;
        k_2 = k*k;

        for (int i = 0; i < MAT_ROWS; i++){
            kxy[i] = ((-(float)MAT_ROWS / 2.0 + i) + 0.5) * delkx; // Center bins
        }

        generate_star_gaussian(software_generated_input_data, MAT_ROWS, sigma, intensity, noise_stddev);

        std::memcpy(bo_input_map, software_generated_input_data, MAT_SIZE * sizeof(cmpx_data_t));
        std::memcpy(bo_kxy_map, kxy, MAT_ROWS * sizeof(float));

        TIMER(cpuTimeStart);
        angular_spectrum_propagation(software_generated_input_data, MAT_ROWS, wavelength, distance, pixel_scale);
        TIMER(cpuTimeEnd);

        cpu_time_arithmetic_mean += ELAPSED(cpuTimeEnd, cpuTimeStart);
        cpu_time_geometric_mean *= ELAPSED(cpuTimeEnd, cpuTimeStart);


        // propagate any argument changes to kerenl
        run.set_arg(1, scale);
        run.set_arg(2, distance);
        run.set_arg(3, k_2);

        TIMER(syncInputTime);
        bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        bo_kxy.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        TIMER(syncInputTimeEnd);

        TIMER(runTime);
        run.start();
        run.wait();
        TIMER(runTimeEnd);

        TIMER(syncOutputTime);
        bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
        TIMER(syncOutputTimeEnd);

        device_exc_time_arithmetic_mean += ELAPSED(runTimeEnd, runTime);
        device_exc_time_geometric_mean  *= ELAPSED(runTimeEnd, runTime);

        device_tot_time_arithmetic_mean += (ELAPSED(syncOutputTimeEnd, syncOutputTime) + ELAPSED(syncInputTimeEnd, syncInputTime) + ELAPSED(runTimeEnd, runTime));
        device_tot_time_geometric_mean *= (ELAPSED(syncOutputTimeEnd, syncOutputTime) + ELAPSED(syncInputTimeEnd, syncInputTime) + ELAPSED(runTimeEnd, runTime));

        printf(timing_summary,  MAT_ROWS, MAT_COLS, 
                                ELAPSED(cpuTimeEnd, cpuTimeStart),
                                ELAPSED(syncOutputTimeEnd, syncOutputTime), 
                                ELAPSED(syncInputTimeEnd, syncInputTime), 
                                ELAPSED(runTimeEnd, runTime),
                                ELAPSED(runTimeEnd, runTime) + ELAPSED(syncOutputTimeEnd, syncOutputTime) + ELAPSED(syncInputTimeEnd, syncInputTime));

        validate_results(software_generated_input_data, bo_output_map, MAT_SIZE, tolerance, false, print_errors);
    }

    device_exc_time_arithmetic_mean /= rounds;
    device_tot_time_arithmetic_mean /= rounds;
    cpu_time_arithmetic_mean /= rounds;

    device_exc_time_geometric_mean = std::pow(device_exc_time_geometric_mean, 1.0 / rounds);
    device_tot_time_geometric_mean = std::pow(device_tot_time_geometric_mean, 1.0 / rounds);;
    cpu_time_geometric_mean = std::pow(cpu_time_geometric_mean, 1.0 / rounds);;

    printf(final_results, MAT_ROWS, MAT_COLS,
                          rounds,

                          device_tot_time_geometric_mean,
                          device_exc_time_geometric_mean,
                          cpu_time_geometric_mean,

                          device_tot_time_arithmetic_mean,
                          device_exc_time_arithmetic_mean,
                          cpu_time_arithmetic_mean);

    return 0;
}
