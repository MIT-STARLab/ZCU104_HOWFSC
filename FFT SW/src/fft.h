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



//////////////////////////////////    1D FFT     //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

/**
 * Implementation of the 1D Cooley-Tukey FFT. Vector Interface.
 * Replaces a by its forward discrete Fourier transform, if invert is input as 0;
 * Replaces a by its inverse discrete Fourier transform times length, if isign is input as 1.
 * Non-recursive (iterative), In-place, O(n log(n)).
 * 
 * @param a  data input vector, must have a size of power of 2.
 * @param invert 0 if forward fft, 1 if inverse fft.
 * @param scale whether to scale results by n if inverse fft.
 */
void fft_vector(vector<cmpx_data_t> &a, bool invert, bool scale = false);


/**
 * Implementation of the 1D Cooley-Tukey FFT. Array Interface.
 * Replaces a by its forward discrete Fourier transform, if invert is input as 0;
 * Replaces a by its inverse discrete Fourier transform times length, if isign is input as 1.
 * Non-recursive (iterative), In-place, O(n log(n)).
 * 
 * @param a  data input array, must have a size of power of 2.
 * @param length of the input array.
 * @param invert 0 if forward fft, 1 if inverse fft.
 * @param scale whether to scale results by n if inverse fft.
 */
void fft(cmpx_data_t *a, int length, bool invert, bool scale = false);


/**
 * Implementation of the 1D Cooley-Tukey FFT. Vector Interface.
 * Replaces a by its forward discrete Fourier transform, if invert is input as 0;
 * Replaces a by its inverse discrete Fourier transform times length, if isign is input as 1.
 * Recursive, O(n log(n)).
 * 
 * @param a  data input vector, must have a size of power of 2.
 * @param invert 0 if forward fft, 1 if inverse fft.
 * @param scale whether to scale results by n if inverse fft.
 */
void fft_recursive(vector<cmpx_data_t> &a, bool invert, bool scale = false);





//////////////////////////////////    2D FFT     //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

/**
 * Implementation of the 2D FFT, uses the 1D implementation of Cooley-Tukey FFT. Array Interface.
 * Replaces a by its forward discrete Fourier transform, if invert is input as 0;
 * Replaces a by its inverse discrete Fourier transform times nm, if isign is input as 1.
 * 
 * Non-recursive (iterative). In-place. O(n^2 log(n))
 * 
 * @param a  2D input array in row-major order.
 * @param n the number of cols in the array, must be a power of 2.
 * @param m the number of rows in the array, must be a power of 2.
 * @param invert 0 if forward fft, 1 if inverse fft.
 * @param scale whether to scale results by nm if inverse fft.
 */
void fft2d(cmpx_data_t *a, int m, int n, bool invert, bool scale = false);




/**
 * Implementation of the 2D FFT, uses the 1D implementation of Cooley-Tukey FFT. Vector Interface.
 * Replaces a by its forward discrete Fourier transform, if invert is input as 0;
 * Replaces a by its inverse discrete Fourier transform times nm, if isign is input as 1.
 * 
 * Non-recursive (iterative). In-place. O(n^2 log(n))
 * 
 * @param a  2D input array, all inner vectors must match in size.
 * @param invert 0 if forward fft, 1 if inverse fft.
 * @param scale whether to scale results by nm if inverse fft.
 */
void fft2d_vector(vector<vector<cmpx_data_t>> &a, bool invert, bool scale = false);




//////////////////////////////    Data Generators    //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

/**
 * Produces FFT Testing Data and expected output
 * 
 * @param input_data   a vector that will later contatins the input data, must be a power of 2 size
 * @param output_data  a vector that will later contatins the results of applying fft on the input data, must match the size of input_data
 * @param frequencies  a vector containing the frequencies to be included
 * @param amplitudes   a vector containing the corresponding amplitudes, must match the size of frequencies
 * @param time_range   it's first element contatins the timing of first measurment. the second the last measurment.
 * @param invert       whether you want the output inversefft or fft
 * @param scale        whether to scale results by n if inverse fft.
 */
void fft_data_generator_vector(vector<cmpx_data_t> &input_data, vector<cmpx_data_t> &output_data, vector<float> &frequencies, vector<float> &amplitudes, vector<float> &time_range, bool invert, bool scale = false);


/**
 * Produces FFT Testing Data and expected output
 * 
 * @param input_data   an array that will later contatins the input data, must be a power of 2 size
 * @param output_data  an array that will later contatins the results of applying fft on the input data, must match the size of input_data
 * @param fft_data_size must be a power of 2.
 * @param frequencies  a vector containing the frequencies to be included
 * @param amplitudes   a vector containing the corresponding amplitudes, must match the size of frequencies
 * @param time_range   it's first element contatins the timing of first measurment. the second the last measurment.
 * @param invert       whether you want the output inversefft or fft
 * @param scale whether to scale results by n if inverse fft.
 */
void fft_data_generator(cmpx_data_t *input_data, cmpx_data_t *output_data, int fft_data_size, vector<float> &frequencies, vector<float> &amplitudes, vector<float> &time_range, bool invert, bool scale = false);


/**
 * Produces random 1D FFT Testing Data ~ Unif(100, 1000) and expected output. Array Interface.
 * @param input_data  input 1d array.
 * @param output_data ouput 1d array.
 * @param m num of data samples, must be a power of 2.
 * @param invert 1 if inverse fft.
 * @param scale whether to scale results by nm if inverse fft.
 */
void fft_random_data_generator(cmpx_data_t *input_data, cmpx_data_t *output_data, int n, bool invert, bool scale = false);


/**
 * Produces random 2D FFT Testing Data ~ Unif(100, 1000) and expected output.
 * The input and output should have dimentions of power of 2.
 * @param input_data  input 2d vector
 * @param output_data output 2d vector, must match the dimentions of the input array
 * @param invert 1 if inverse 2dfft
 * @param scale whether to scale results by nm if inverse fft.
 */
void fft2d_random_data_generator_vector(vector<vector<cmpx_data_t>> &input_data, vector<vector<cmpx_data_t>> &output_data, bool invert, bool scale = false);


/**
 * Produces random 2D FFT Testing Data ~ Unif(100, 1000) and expected output. Array Interface.
 * @param input_data  input 2d array in row-major order.
 * @param output_data ouput 2d array in row-major order.
 * @param m num of rows, must be a power of 2.
 * @param n num of cols, must be a power of 2.
 * @param invert 1 if inverse 2dfft
 * @param scale whether to scale results by nm if inverse fft.
 */
void fft2d_random_data_generator(cmpx_data_t *input_data, cmpx_data_t *output_data, int m, int n, bool invert, bool scale = false);