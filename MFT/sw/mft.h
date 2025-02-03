/*
 * MIT STAR Lab
 * M.Subhi Abo Rdan (msubhi_a@mit.edu)
 * Last modified on Feb 3, 2024
 * Pure software implementations of Matrix Fourier Transform to help test hardware kernels
 */

#ifndef MFT_SW_H
#define MFT_SW_H

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

/**
 * @param wavefront             complex 2D array representing the pupil plane wavefront
 * @param psf_pixelscale_lamD   the pixelscale of the desired focal plane wavefront in terms of lambda/D
 * @param npsf                  the size of the desired focal plane in pixels
 * 
 * On return, wavefront has the complex wavefront at the focal plane
 */
template <typename T>
void mft(std::complex<T>* wavefront, int npix, int npsf, T psf_pixelscale_lamD, bool forward);


template <typename T>
void complex_dmatmul(std::complex<T> *C, std::complex<T> *A, std::complex<T> *B, int M, int N, int P);


#endif // MFT_SW_H