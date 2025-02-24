/*
 * MIT STAR Lab
 * Pure software (C) implementations of Matrix Fourier Transform to help test hardware kernels
 * 
 * Nicholas Belsten.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>

void mft_forward(int npix, double* *wavefront, double psf_pixelscale_lamD, int npsf);
void mft_reverse(int npsf, double* *wavefront, double psf_pixelscale_lamD, int npix);
void sequential_complex_dmatmul(double *C, double *A, double *B, int M, int N, int P);
