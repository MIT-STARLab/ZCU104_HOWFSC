/*
 *
 * MIT STAR Lab
 * Modified by Subhi in Jun 13
 * FFT Implenetation in C++ to test our Hardware Implementation
 * 
 */

#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <fstream>



using namespace std;


typedef complex<float> cmpx_data_t;


/**
 * Implementation of the 1D Cooley-Tukey FFT, the input should have a length of power of 2. 
 * @param a  input vector
 * @param invert 1 if inverse fft
 */
void fft(vector<cmpx_data_t> &a, bool invert);



/**
 * Produces FFT Testing Data and expected output
 * @param input_data   a vector that will later contatins the input data, must be a power of 2 size
 * @param output_data  a vector that will later contatins the results of applying fft on the input data, must match the size of input_data
 * @param frequencies  a vector containing the frequencies to be included
 * @param amplitudes   a vector containing the corresponding amplitudes, must match the size of frequencies
 * @param time_range   it's first element contatins the timing of first measurment. the second the last measurment.
 * @param invert       whether you want the output inversefft or fft
 */
void fft_data_generator(vector<cmpx_data_t> &input_data, vector<cmpx_data_t> &output_data, vector<float> &frequencies, vector<float> &amplitudes, vector<float> &time_range, bool invert);


