/*
 *
 * MIT STAR Lab
 * Modified by Subhi in Jun 13
 * Pure software implementations of useful FFT function to help test hardware kernels
 * 
 */


#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>

using namespace std;

typedef complex<float> cmpx_data_t;


/**
 * Implementation of the 1D Cooley-Tukey FFT. 
 * Non-recursive (iterative).
 * In-place.
 * O(n log(n))
 * 
 * @param a  data input vector, must have a size of power of 2.
 * @param invert 0 if forward fft, 1 if inverse fft.
 */
void fft(vector<cmpx_data_t> &a, bool invert);


/**
 * Implementation of the 1D Cooley-Tukey FFT. 
 * Recursive.
 * In-place.
 * O(n log(n))
 * 
 * @param a  data input vector, must have a size of power of 2.
 * @param invert 0 if forward fft, 1 if inverse fft.
 */
void fft_recursive(vector<cmpx_data_t> &a, bool invert);


/**
 * Implementation of the 2D FFT, uses the 1D implementation of Cooley-Tukey FFT.
 * The input must have dimentions of power of 2.
 * Non-recursive (iterative).
 * In-place.
 * O(n^2 log(n))
 * 
 * @param a  2D input array, all inner vectors must match in size.
 * @param invert 0 if forward fft, 1 if inverse fft.
 */
void fft2d(vector<vector<cmpx_data_t>> &a, bool invert);


/**
 * Produces FFT Testing Data and expected output
 * 
 * @param input_data   a vector that will later contatins the input data, must be a power of 2 size
 * @param output_data  a vector that will later contatins the results of applying fft on the input data, must match the size of input_data
 * @param frequencies  a vector containing the frequencies to be included
 * @param amplitudes   a vector containing the corresponding amplitudes, must match the size of frequencies
 * @param time_range   it's first element contatins the timing of first measurment. the second the last measurment.
 * @param invert       whether you want the output inversefft or fft
 */
void fft_data_generator(vector<cmpx_data_t> &input_data, vector<cmpx_data_t> &output_data, vector<float> &frequencies, vector<float> &amplitudes, vector<float> &time_range, bool invert);



/**
 * Produces random FFT Testing Data and expected outpu
 * The input and output should have dimentions of power of 2.
 * @param input_data  input array
 * @param output_data input array, must match the dimentions of the input array
 * @param invert 1 if inverse fft
 */
void fft2d_random_data_generator(vector<vector<cmpx_data_t>> &input_data, vector<vector<cmpx_data_t>> &output_data, bool invert);

/**
 * Same as above but an array interface of the matrix with row-major order
 */
void fft2d_random_data_generator_array(cmpx_data_t *input_data, cmpx_data_t *output_data, int m, int n, bool invert);
