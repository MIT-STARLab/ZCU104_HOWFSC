/*
 * MIT STAR Lab
 * M.Subhi Abo Rdan (msubhi_a@mit.edu)
 * Last modified in Dec 29, 2024
 */

#ifndef ANG_SPEC_PROP_H
#define ANG_SPEC_PROP_H

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include "fft.h"

/**
 * @brief Performs angular spectrum propagation of a 2D wavefront.
 * 
 * This function propagates a 2D wavefront in free space over a specified distance using the angular 
 * spectrum method. It performs forward and inverse Fourier transforms to switch between the spatial 
 * and frequency domains, applying a transfer function to model free-space propagation.
 * 
 * @param wavefront     Input/output array of size n * n in row-major order. On input, it represents 
 *                      the complex-valued 2D wavefront in the spatial domain. On output, it contains 
 *                      the propagated wavefront in the spatial domain after propagation.
 * @param n             The size of the square wavefront grid (number of rows and columns). Must be 
 *                      a power of 2 for FFT compatibility.
 * @param wavelength    The wavelength of the light (or wave) in meters.
 * @param distance      The propagation distance in meters.
 * @param pixel_scale   The physical size of each pixel in the grid in units of [pixels/meter].
 *                      Determines the spatial resolution of the input wavefront.
 */
extern "C" void angular_spectrum_propagation(cmpx_data_t* wavefront, int n, float wavelength, float distance, float pixel_scale);


// Helpers
extern "C" void construct_ishifted_transform_function(cmpx_data_t* tf, int n, float wavelength, float distance, float pixel_scale);
extern "C" void construct_transform_function(cmpx_data_t* tf,int n, float wavelength, float distance, float pixel_scale);
extern "C" void elementwise_cmpx_mul(cmpx_data_t* arr1, cmpx_data_t* arr2, int m, int n);


// Testing Data Helper
extern "C" void generate_star_gaussian(cmpx_data_t* arr, int size, double sigma, double intensity, double noise_stddev, bool noise = true);

#endif