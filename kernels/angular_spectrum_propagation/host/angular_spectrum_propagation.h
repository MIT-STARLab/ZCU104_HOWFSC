/*
 * MIT STAR Lab
 * M.Subhi Abo Rdan (msubhi_a@mit.edu)
 * Last modified on January 15, 2025
 * Pure software implementations of Angular Spectrum Propagation Method to help test hardware kernels
 */

#ifndef ANG_SPEC_PROP_SW_H
#define ANG_SPEC_PROP_SW_H

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include "fft.h"

/**
 * @brief Performs angular spectrum propagation on a 2D wavefront.
 * 
 * This function propagates a 2D wavefront in free space over a specified distance using the angular 
 * spectrum method. It utilizes forward and inverse Fourier transforms to switch between the spatial 
 * and frequency domains and applies a transfer function for free-space propagation.
 * 
 * @param wavefront     Input/output array of size n * n in row-major order. On input, it represents 
 *                      the complex-valued 2D wavefront in the spatial domain. On output, it contains 
 *                      the propagated wavefront after propagation.
 * @param n             The size of the square wavefront grid (number of rows and columns). Must be 
 *                      a power of 2 for FFT compatibility.
 * @param wavelength    The wavelength of the wave in meters.
 * @param distance      The propagation distance in meters.
 * @param pixel_scale   The physical size of each pixel in the grid in units of [pixels/meter].
 *                      Determines the spatial resolution of the wavefront.
 */
template <typename T>
void angular_spectrum_propagation(std::complex<T>* wavefront, int n, T wavelength, T distance, T pixel_scale);


/**
 * @brief Propagates the wavefront using a shifted method.
 */
template <typename T>
void propagate_ishifted(std::complex<T>* wavefront, int n, T wavelength, T distance, T pixel_scale, bool mod = false);

/**
 * @brief Propagates the wavefront using a basic method.
 */
template <typename T>
void propagate(std::complex<T>* wavefront, int n, T wavelength, T distance, T pixel_scale, bool mod = false);

/**
 * @brief Constructs the transfer function for propagation.
 */
template <typename T>
void construct_transform_function(std::complex<T>* tf, int n, T wavelength, T distance, T pixel_scale, bool mod = false);

/**
 * @brief Constructs a shifted transfer function for propagation.
 */
template <typename T>
void construct_ishifted_transform_function(std::complex<T>* tf, int n, T wavelength, T distance, T pixel_scale, bool mod = false);

/**
 * @brief Generates a Gaussian wavefront with noise.
 */
template <typename T>
void generate_star_gaussian(std::complex<T>* arr, int size, T sigma, T intensity, T noise_stddev, bool noise = true);


#ifdef PYTHON_TESTING_FLOATS
extern "C" {
    void angular_spectrum_propagation_float_wrapper(std::complex<float>* wavefront, int n, float wavelength, float distance, float pixel_scale) {
        angular_spectrum_propagation(wavefront, n, wavelength, distance, pixel_scale);
    }
    void propagate_ishifted_float_wrapper(std::complex<float>* wavefront, int n, float wavelength, float distance, float pixel_scale, bool mod = false){
        propagate_ishifted(wavefront, n, wavelength, distance, pixel_scale, mod);
    }
    void propagate_float_wrapper(std::complex<float>* wavefront, int n, float wavelength, float distance, float pixel_scale, bool mod = false){
        propagate(wavefront, n, wavelength, distance, pixel_scale, mod);
    }
    void construct_transform_function_float_wrapper(std::complex<float>* tf, int n, float wavelength, float distance, float pixel_scale, bool mod = false){
        construct_transform_function(tf, n, wavelength, distance, pixel_scale, mod);
    }
    void construct_ishifted_transform_function_float_wrapper(std::complex<float>* tf, int n, float wavelength, float distance, float pixel_scale, bool mod = false){
        construct_ishifted_transform_function(tf, n, wavelength, distance, pixel_scale, mod);
    }
    void generate_star_gaussian_float_wrapper(std::complex<float>* arr, int size, float sigma, float intensity, float noise_stddev, bool noise = true){
        generate_star_gaussian(arr, size, sigma, intensity, noise_stddev, noise);
    }
}
#endif // PYTHON_TESTING_FLOATS


#ifdef PYTHON_TESTING_DOUBLES
extern "C" {
    void angular_spectrum_propagation_double_wrapper(std::complex<double>* wavefront, int n, double wavelength, double distance, double pixel_scale) {
        angular_spectrum_propagation(wavefront, n, wavelength, distance, pixel_scale);
    }
    void propagate_ishifted_double_wrapper(std::complex<double>* wavefront, int n, double wavelength, double distance, double pixel_scale, bool mod = false){
        propagate_ishifted(wavefront, n, wavelength, distance, pixel_scale, mod);
    }
    void propagate_double_wrapper(std::complex<double>* wavefront, int n, double wavelength, double distance, double pixel_scale, bool mod = false){
        propagate(wavefront, n, wavelength, distance, pixel_scale, mod);
    }
    void construct_transform_function_double_wrapper(std::complex<double>* tf, int n, double wavelength, double distance, double pixel_scale, bool mod = false){
        construct_transform_function(tf, n, wavelength, distance, pixel_scale, mod);
    }
    void construct_ishifted_transform_function_double_wrapper(std::complex<double>* tf, int n, double wavelength, double distance, double pixel_scale, bool mod = false){
        construct_ishifted_transform_function(tf, n, wavelength, distance, pixel_scale, mod);
    }
    void generate_star_gaussian_double_wrapper(std::complex<double>* arr, int size, double sigma, double intensity, double noise_stddev, bool noise = true){
        generate_star_gaussian(arr, size, sigma, intensity, noise_stddev, noise);
    }
}
#endif // PYTHON_TESTING_DOUBLES


#endif // ANG_SPEC_PROP_SW_H