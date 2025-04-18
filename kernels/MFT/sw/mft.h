/*
 * MIT STAR Lab
 * Pure software implementations of Matrix Fourier Transform to help test hardware kernels (C++)
 * 
 * C code inherted from Nicholas Belsten.
 * Last modified on Feb 23, 2025 
 * by M.Subhi Abo Rdan (msubhi_a@mit.edu)
 */

#ifndef MFT_SW_H
#define MFT_SW_H

#include <complex>

/**
 * @brief Computes the forward modal Fourier transform (MFT).
 *
 * @tparam T Data type of the complex matrix elements (e.g., float, double).
 * @param wavefront The pupil plane wavefront, stored in row-major order.
 * @param psf_pixelscale_lamD The pixel scale of the desired focal plane wavefront in terms of λ/D.
 * @param npsf The size of the desired focal plane in pixels.
 * @param npix The size of the input pupil array.
 *
 * Returns
 * -------
 * The complex wavefront at the focal plane, stored in row-major order.
 *
 * @note The `wavefront` array is deallocated and reallocated with a new size.
 */
template <typename T>
void mft_forward(std::complex<T> **wavefront, T psf_pixelscale_lamD, int npsf, int npix);
template <typename T>
void mft_reverse(std::complex<T> **wavefront, T psf_pixelscale_lamD, int npsf, int npix);

/**
 * @brief Performs matrix multiplication for complex matrices using recursive blocking.
 *
 * This function multiplies two complex matrices, A (M×N) and B (N×P), 
 * and stores the result in C (M×P). It leverages recursive partitioning 
 * and OpenMP parallelization for improved performance.
 *
 * @tparam T Data type of the complex matrix elements (e.g., float, double).
 * @param C Pointer to the output matrix (M×P), stored in a contiguous 1D array.
 * @param A Pointer to the first input matrix (M×N), stored in a contiguous 1D array.
 * @param B Pointer to the second input matrix (N×P), stored in a contiguous 1D array.
 * @param M Number of rows in A and C.
 * @param N Number of columns in A and rows in B.
 * @param P Number of columns in B and C.
 */
template <typename T>
void complex_mat_mul(std::complex<T> *C, std::complex<T> *A, std::complex<T> *B, int M, int N, int P);

#endif // MFT_SW_H