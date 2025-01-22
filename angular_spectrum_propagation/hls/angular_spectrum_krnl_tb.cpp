/*
 * MIT STAR Lab
 * M.Subhi Abo Rdan (msubhi_a@mit.edu)
 * Last modified on January 16, 2024
 * Pure software implementations of Angular Spectrum Propagation Method to help test hardware kernels
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>

#include "angular_spectrum_krnl.h"
#include "angular_spectrum_propagation.h"

using namespace std;

/**
 * Checks if two complex numbers are equal within a tolarance
 */
bool approx_equal(const cmpx_data_t &a, const cmpx_data_t &b, double tolerance) {

    bool real_compare = (a.real()<tolerance && b.real()<tolerance) || (std::abs((a.real() - b.real())/std::max(std::abs(a.real()), std::abs(b.real()))) < tolerance);
    bool imag_compare = (a.imag()<tolerance && b.imag()<tolerance) || (std::abs((a.imag() - b.imag())/std::max(std::abs(a.imag()), std::abs(b.imag()))) < tolerance);

    return  real_compare && imag_compare;
}

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


/**
 * Validate results
 */
void validate_results(cmpx_data_t *expected_output, cmpx_data_t *real_output, int data_size, double tolerance){

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

}


cmpx_data_t software_generated_data[MAT_SIZE];

cmpx_data_t hardware_input_data[MAT_SIZE];
cmpx_data_t hardware_output_data[MAT_SIZE];
data_t kxy[MAT_ROWS];

int main() {
    cout << "INFO:  starting angular spectrum propagation kernal testbench" << endl;

    cout << "INFO:  Generate Testing Data using Pure software implementaion" << endl;
    //// Testing parameters
    data_t sigma        = 10;
    data_t intensity    = 1e5;
    data_t noise_stddev = 0.01;
    bool   noise        = true;

    data_t wavelength   = 500e-9;
    data_t distance     = 1000e-3;
    data_t pixel_scale  = 10e-3 / ((data_t)MAT_ROWS / 2);

    data_t delkx = 2.0 * M_PI / (pixel_scale * MAT_ROWS);
    data_t k = 2.0 * M_PI / wavelength;
    data_t k_2 = k*k;
    data_t scale = (1/(data_t) MAT_SIZE);

    generate_star_gaussian(software_generated_data, MAT_ROWS, sigma, intensity, noise_stddev, noise);
    cout << STR_PASSED <<"generate_star_gaussian" << endl;
    std::memcpy(hardware_input_data, software_generated_data, sizeof(hardware_input_data));
    cout << STR_PASSED <<"memcpy" << endl;
    angular_spectrum_propagation(software_generated_data, MAT_ROWS, wavelength, distance, pixel_scale);
    cout << "Testing data generated" << endl;


    cout << "Call the angular spectrum propagation Kernel" << endl;
    for (int i = 0; i < MAT_ROWS; i++){
        kxy[i] = ((-(data_t)MAT_ROWS / 2.0 + i) + 0.5) * delkx; // Center bins
    }

    angular_spectrum( // krnl
        true,   // direction
        scale,
        distance,
        k_2,
        kxy,
        hardware_input_data,
        hardware_output_data
        );

    cout << "Saving Results to files <software_output.txt> and <hardware_output.txt>" << endl;
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


    cout << "Validating Results" << endl;
    validate_results(software_generated_data, hardware_output_data, MAT_SIZE, 1e-2);

    return 0;
}