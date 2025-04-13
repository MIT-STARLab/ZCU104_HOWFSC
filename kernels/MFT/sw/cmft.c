/*
 * MIT STAR Lab
 * Pure software (C) implementations of Matrix Fourier Transform to help test hardware kernels
 * 
 * Nicholas Belsten.
 */

#include "cmft.h"

#define NR_END 1
#define FREE_ARG char*


void nrerror(const char *message) {
    fprintf(stderr, "Error: %s\n", message);
    exit(EXIT_FAILURE);
}


void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}


void sequential_complex_dmatmul(double *C, double *A, double *B, int M, int N, int P) {
    int i, j, k;
    double a_real, a_im, b_real, b_im, c_real, c_im;

    // Each element in A, B, C consists of 2 doubles (real and imaginary parts)
    for (i = 0; i < M; i++) {
        for (j = 0; j < P; j++) {
            double real_sum = 0.0;
            double im_sum = 0.0;
            for (k = 0; k < N; k++) {
                // Real and imaginary parts of A[i][k]
                a_real = A[2 * (i * N + k) + 1];     // 2 * (i * N + k) points to the real part
                a_im = A[2 * (i * N + k) + 2];   // 2 * (i * N + k) + 1 points to the imaginary part

                // Real and imaginary parts of B[k][j]
                b_real = B[2 * (k * P + j) + 1];     // 2 * (k * P + j)
                b_im = B[2 * (k * P + j) + 2];

                // Complex multiplication: (a_real + i*a_im) * (b_real + i*b_im)
                real_sum += a_real * b_real - a_im * b_im;
                im_sum += a_real * b_im + a_im * b_real;
            }
            // Store the result in C[i][j]
            C[2 * (i * P + j) + 1] = real_sum;    // Real part
            C[2 * (i * P + j) + 2] = im_sum;  // Imaginary part
        }
    }
}


// Parameters
// ----------
// pupil : complex 2D array
// the pupil plane wavefront
// psf_pixelscale_lamD : scalar
// the pixelscale of the desired focal plane wavefront in terms
// of lambda/D
// npsf : integer
// the size of the desired focal plane in pixels

// Returns
// -------
// complex 2D array
// the complex wavefront at the focal plane
// 
// Note: The wavefront is deallocated then reallocated at a different size in this function
// """

void mft_forward(int npix, double* *wavefront, double psf_pixelscale_lamD, int npsf){
    double dx = 1.0 / (double)npix;
    double du = psf_pixelscale_lamD;

    //Xs = (xp.arange(npix, dtype=float) - (npix / 2)) * dx
    double Xs[npix];
    for (int i = 0; i < npix; i++) {
        //Xs[i] = ((-(double)npix / 2.0 + i) + 0.5) * dx; //Note to Kian, off by 1/2?
        Xs[i] = (((double)i - (double)npix/2.0)) * dx;
        //printf("%f\n",Xs[i]);
    }

    // Us = (xp.arange(npsf, dtype=float) - npsf / 2) * du
    double Us[npsf];
    for (int i = 0; i < npsf; i++) {
        //Us[i] = ((-(double)npsf / 2.0 + i) + 0.5) * du; //Note to Kian, off by 1/2?
        Us[i] = (((double)i - (double)npsf/2.0)) * du;
        //printf("%f\n",Us[i]);
    }

    double* Mx = dvector(1,2*npsf*npix); //Needs freeing. is freed at end of function. Shape is npsf by npix
    double* My = dvector(1,2*npsf*npix); //Needs freeing. is freed at end of function. Shape is npix by npsf

    //Outer of Us, Xs
    //Result shoudl be size [npsf,npix]
    for(int i = 0; i < npsf; i++){
        for(int j = 0; j < npix; j++){ //Iterate first across rows for row major order
            //We are at index [i, j] in zero based
            int index = (i*npix + j)*2 + 1; //1 based index, this is the real part. Note this index is just for storing result
            double x_outer = Us[i] * Xs[j];
            double complex x = cexp(-2*M_PI*I * x_outer) / npix;
            Mx[index] = creal(x);
            Mx[index+1] = cimag(x);
        }
    }

    //Outer of Xs, Us
    //Result should be size [npix, npsf]
    for(int i = 0; i < npix; i++){
        for(int j = 0; j < npsf; j++){ //Iterate first across rows for row major order
            //We are at index [i, j] in zero based
            int index = (i*npsf + j)*2 + 1; //1 based index, this is the real part
            double y_outer = Xs[i] * Us[j];
            double complex y = cexp(-2*M_PI*I * y_outer) * du;
            My[index] = creal(y);
            My[index+1] = cimag(y);
        }
    }

    //Need to do complex matmul here
    //Mx@pupil@My
    double* temp = dvector(1,2*npix*npsf); //is freed at bottom of function

    sequential_complex_dmatmul(temp,*wavefront,My, npix, npix, npsf); // Result is size M by P = npix by npsf

    free_dvector(*wavefront,1,npix*npix*2);
    *wavefront = dvector(1,npsf*npsf*2);

    sequential_complex_dmatmul(*wavefront,Mx,temp, npsf, npix, npsf); // B is npix by npsf = N by P. Result is M by P is npsf by npsf
    free_dvector(temp,1,2*npix*npsf);
    free_dvector(Mx,1,2*npsf*npix);
    free_dvector(My,1,2*npsf*npix);
}




// Note: The wavefront is deallocated then reallocated at a different size in this function
void mft_reverse(int npsf, double* *wavefront, double psf_pixelscale_lamD, int npix){
    double dx = 1.0 / (double)npix;
    double du = psf_pixelscale_lamD;

    double Xs[npix];
    for (int i = 0; i < npix; i++) {
        //Xs[i] = ((-(double)npix / 2.0 + i) + 0.5) * dx; //Note to Kian, off by 1/2?
        Xs[i] = (((double)i - (double)npix/2.0)) * dx;
        //printf("%f\n",Xs[i]);
    }

    double Us[npsf];
    for (int i = 0; i < npsf; i++) {
        //Us[i] = ((-(double)npsf / 2.0 + i) + 0.5) * du; //Note to Kian, off by 1/2?
        Us[i] = (((double)i - (double)npsf/2.0)) * du;
        //printf("%f\n",Us[i]);
    }

    double* Mx = dvector(1,2*npsf*npix); //Needs freeing. Shape is npsf by npix
    double* My = dvector(1,2*npsf*npix); //Needs freeing. Shape is npix by npsf

    // ux = xp.outer(Xs, Us)
    // yv = xp.outer(Us, Xs)

    // My = xp.exp(-1j*2*np.pi*yv) 
    // Mx = xp.exp(-1j*2*np.pi*ux) 

    //Outer of Us, Xs
    //Result shoudl be size [npsf,npix]
    for(int i = 0; i < npsf; i++){
        for(int j = 0; j < npix; j++){ //Iterate first across rows for row major order
            //We are at index [i, j] in zero based
            int index = (i*npix + j)*2 + 1; //1 based index, this is the real part. Note this index is just for storing result
            double yv = Us[i] * Xs[j];
            double complex y = cexp(-2*M_PI*I * yv) / npix;
            My[index] = creal(y);
            My[index+1] = cimag(y);
        }
    }

    //Outer of Xs, Us
    //Result should be size [npix, npsf]
    for(int i = 0; i < npix; i++){
        for(int j = 0; j < npsf; j++){ //Iterate first across rows for row major order
            //We are at index [i, j] in zero based
            int index = (i*npsf + j)*2 + 1; //1 based index, this is the real part
            double ux = Xs[i] * Us[j];
            double complex x = cexp(-2*M_PI*I * ux) * du;
            Mx[index] = creal(x);
            Mx[index+1] = cimag(x);
        }
    }

    //return Mx@fpwf@My
    double* temp = dvector(1,2*npix*npsf);

    sequential_complex_dmatmul(temp,*wavefront,My, npsf, npsf, npix); //B is [npsf,npix], A is [npsf,npsf], C is [npsf,npix]
    sequential_complex_dmatmul(*wavefront,Mx,temp, npix, npsf, npix); //B is [npsf,npix], A is [npix,npsf], C is [npix,npix]

    free_dvector(Mx,1,2*npsf*npix);
    free_dvector(My,1,2*npsf*npix);
    free_dvector(temp,1,2*npix*npsf);
}

