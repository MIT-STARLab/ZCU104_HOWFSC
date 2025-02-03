/*
 * MIT STAR Lab
 * M.Subhi Abo Rdan (msubhi_a@mit.edu)
 * Last modified on Feb 3, 2024
 * Pure software implementations of Matrix Fourier Transform to help test hardware kernels
 * Code inherted from Nick.
 */

#include "mft.h"


template <typename T>
void mft(std::complex<T>* wavefront, int npix, int npsf, T psf_pixelscale_lamD, bool forward)
{
    T dx = 1.0 / (T)npix;
    T du = psf_pixelscale_lamD;

    //Xs = (xp.arange(npix, dtype=float) - (npix / 2)) * dx
    T Xs[npix];
    for (int i = 0; i < npix; i++) {
        Xs[i] = (((T)i - (T)npix/2.0)) * dx;
    }

    // Us = (xp.arange(npsf, dtype=float) - npsf / 2) * du
    T Us[npsf];
    for (int i = 0; i < npsf; i++) {
        Us[i] = (((T)i - (T)npsf/2.0)) * du;
    }

    std::complex<T>* Mx = (std::complex<T>*) malloc(npsf * npix * sizeof(std::complex<T>)); //Shape is npsf by npix
    std::complex<T>* My = (std::complex<T>*) malloc(npsf * npix * sizeof(std::complex<T>)); //Shape is npix by npsf

    //Outer of Us, Xs
    //Result shoudl be size [npsf,npix]
    for(int i = 0; i < npsf; i++){
        for(int j = 0; j < npix; j++){  //Iterate first across rows for row major order
            //We are at index [i, j] in zero based
            int index = (i*npix + j) + 1; //1 based index, this is the real part. Note this index is just for storing result
            T x_outer = Us[i] * Xs[j];
            complex<T> x = cexp(-2*M_PI*I * x_outer) / npix;
            Mx[index] = x;
        }
    }

    //Outer of Xs, Us
    //Result should be size [npix, npsf]
    for(int i = 0; i < npix; i++){
        for(int j = 0; j < npsf; j++){ //Iterate first across rows for row major order
            //We are at index [i, j] in zero based
            int index = (i*npsf + j) + 1; //1 based index, this is the real part
            T y_outer = Xs[i] * Us[j];
            complex<T> y = cexp(-2*M_PI*I * y_outer) * du;
            My[index] = y;
        }
    }

    //Need to do complex matmul here
    //Mx@pupil@My
    std::complex<T>* temp = (std::complex<T>*) malloc(npsf * npix * sizeof(std::complex<T>));


    // complex_dmatmul(temp,*wavefront,My, npix, npix, npsf); // Result is size M by P = npix by npsf

    // free_dvector(*wavefront,1,npix*npix*2);
    // *wavefront = dvector(1,npsf*npsf*2);

    // complex_dmatmul(*wavefront,Mx,temp, npsf, npix, npsf); // B is npix by npsf = N by P. Result is M by P is npsf by npsf


    free(temp);
    free(Mx);
    free(My);
}





// // Parameters
// // ----------
// // pupil : complex 2D array
// // the pupil plane wavefront
// // psf_pixelscale_lamD : scalar
// // the pixelscale of the desired focal plane wavefront in terms
// // of lambda/D
// // npsf : integer
// // the size of the desired focal plane in pixels

// // Returns
// // -------
// // complex 2D array
// // the complex wavefront at the focal plane
// // 
// // Note: The wavefront is deallocated then reallocated at a different size in this function
// // """
// void mft_forward(int npix, double* *wavefront, double psf_pixelscale_lamD, int npsf){
//     double dx = 1.0 / (double)npix;
//     double du = psf_pixelscale_lamD;

//     //Xs = (xp.arange(npix, dtype=float) - (npix / 2)) * dx
//     double Xs[npix];
//     for (int i = 0; i < npix; i++) {
//         //Xs[i] = ((-(double)npix / 2.0 + i) + 0.5) * dx; //Note to Kian, off by 1/2?
//         Xs[i] = (((double)i - (double)npix/2.0)) * dx;
//         //printf("%f\n",Xs[i]);
//     }

//     // Us = (xp.arange(npsf, dtype=float) - npsf / 2) * du
//     double Us[npsf];
//     for (int i = 0; i < npsf; i++) {
//         //Us[i] = ((-(double)npsf / 2.0 + i) + 0.5) * du; //Note to Kian, off by 1/2?
//         Us[i] = (((double)i - (double)npsf/2.0)) * du;
//         //printf("%f\n",Us[i]);
//     }

//     double* Mx = dvector(1,2*npsf*npix); //Needs freeing. is freed at end of function. Shape is npsf by npix
//     double* My = dvector(1,2*npsf*npix); //Needs freeing. is freed at end of function. Shape is npix by npsf

//     //Outer of Us, Xs
//     //Result shoudl be size [npsf,npix]
//     for(int i = 0; i < npsf; i++){
//         for(int j = 0; j < npix; j++){ //Iterate first across rows for row major order
//             //We are at index [i, j] in zero based
//             int index = (i*npix + j)*2 + 1; //1 based index, this is the real part. Note this index is just for storing result
//             double x_outer = Us[i] * Xs[j];
//             complex x = cexp(-2*M_PI*I * x_outer) / npix;
//             Mx[index] = creal(x);
//             Mx[index+1] = cimag(x);
//         }
//     }

//     //Outer of Xs, Us
//     //Result should be size [npix, npsf]
//     for(int i = 0; i < npix; i++){
//         for(int j = 0; j < npsf; j++){ //Iterate first across rows for row major order
//             //We are at index [i, j] in zero based
//             int index = (i*npsf + j)*2 + 1; //1 based index, this is the real part
//             double y_outer = Xs[i] * Us[j];
//             complex y = cexp(-2*M_PI*I * y_outer) * du;
//             My[index] = creal(y);
//             My[index+1] = cimag(y);
//         }
//     }

//     //Need to do complex matmul here
//     //Mx@pupil@My
//     double* temp = dvector(1,2*npix*npsf); //is freed at bottom of function

//     complex_dmatmul(temp,*wavefront,My, npix, npix, npsf); // Result is size M by P = npix by npsf

//     free_dvector(*wavefront,1,npix*npix*2);
//     *wavefront = dvector(1,npsf*npsf*2);

//     complex_dmatmul(*wavefront,Mx,temp, npsf, npix, npsf); // B is npix by npsf = N by P. Result is M by P is npsf by npsf
//     free_dvector(temp,1,2*npix*npsf);
//     free_dvector(Mx,1,2*npsf*npix);
//     free_dvector(My,1,2*npsf*npix);
// }





// // Note: The wavefront is deallocated then reallocated at a different size in this function
// void mft_reverse(int npsf, double* *wavefront, double psf_pixelscale_lamD, int npix){
//     double dx = 1.0 / (double)npix;
//     double du = psf_pixelscale_lamD;

//     double Xs[npix];
//     for (int i = 0; i < npix; i++) {
//         //Xs[i] = ((-(double)npix / 2.0 + i) + 0.5) * dx; //Note to Kian, off by 1/2?
//         Xs[i] = (((double)i - (double)npix/2.0)) * dx;
//         //printf("%f\n",Xs[i]);
//     }

//     double Us[npsf];
//     for (int i = 0; i < npsf; i++) {
//         //Us[i] = ((-(double)npsf / 2.0 + i) + 0.5) * du; //Note to Kian, off by 1/2?
//         Us[i] = (((double)i - (double)npsf/2.0)) * du;
//         //printf("%f\n",Us[i]);
//     }

//     double* Mx = dvector(1,2*npsf*npix); //Needs freeing. Shape is npsf by npix
//     double* My = dvector(1,2*npsf*npix); //Needs freeing. Shape is npix by npsf

//     // ux = xp.outer(Xs, Us)
//     // yv = xp.outer(Us, Xs)

//     // My = xp.exp(-1j*2*np.pi*yv) 
//     // Mx = xp.exp(-1j*2*np.pi*ux) 

//     //Outer of Us, Xs
//     //Result shoudl be size [npsf,npix]
//     for(int i = 0; i < npsf; i++){
//         for(int j = 0; j < npix; j++){ //Iterate first across rows for row major order
//             //We are at index [i, j] in zero based
//             int index = (i*npix + j)*2 + 1; //1 based index, this is the real part. Note this index is just for storing result
//             double yv = Us[i] * Xs[j];
//             complex y = cexp(-2*M_PI*I * yv) / npix;
//             My[index] = creal(y);
//             My[index+1] = cimag(y);
//         }
//     }

//     //Outer of Xs, Us
//     //Result should be size [npix, npsf]
//     for(int i = 0; i < npix; i++){
//         for(int j = 0; j < npsf; j++){ //Iterate first across rows for row major order
//             //We are at index [i, j] in zero based
//             int index = (i*npsf + j)*2 + 1; //1 based index, this is the real part
//             double ux = Xs[i] * Us[j];
//             complex x = cexp(-2*M_PI*I * ux) * du;
//             Mx[index] = creal(x);
//             Mx[index+1] = cimag(x);
//         }
//     }

//     //return Mx@fpwf@My
//     double* temp = dvector(1,2*npix*npsf);

//     complex_dmatmul(temp,*wavefront,My, npsf, npsf, npix); //B is [npsf,npix], A is [npsf,npsf], C is [npsf,npix]
//     complex_dmatmul(*wavefront,Mx,temp, npix, npsf, npix); //B is [npsf,npix], A is [npix,npsf], C is [npix,npix]

//     free_dvector(Mx,1,2*npsf*npix);
//     free_dvector(My,1,2*npsf*npix);
//     free_dvector(temp,1,2*npix*npsf);
// }



template <typename T>
void base_complex_dmatmul(std::complex<T>  **C, std::complex<T>  **A, std::complex<T>  **B, int M, int N, int P) 
{
    int i, j, k;
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= P; j++) {
            std::complex<T> sum(0.0, 0.0);
            for (k = 1; k <= N; k++) {
                sum += (A[i][k] * B[k][j]);
            }
            C[i][j] += sum;
        }
    }
}


// Helper to create a submatrix pointer array for complex matrices
template <typename T>
std::complex<T> **create_complex_submatrix(std::complex<T> **matrix, int row_offset, int col_offset, int rows, int cols) 
{
    std::complex<T> **submatrix = (std::complex<T> **)malloc((rows + 1) * sizeof(std::complex<T> *));
    if (submatrix == NULL) {
        perror("Failed to allocate submatrix pointers");
        exit(EXIT_FAILURE);
    }

    // Each row starts at 2 * ((i - 1) * cols + col_offset) + 1 to account for complex layout
    for (int i = 1; i <= rows; i++) {
        submatrix[i] = &matrix[row_offset + i][col_offset*2]; //plus 1 here? no I don't think so because we'll be calling submatrix[i][1], so we want this to reference submatrix[i][0]

                // Debugging: Print the base address of each row in the submatrix
        // printf("Submatrix row %d points to matrix[%d][%d] -> Real part address: %p\n",
        //         i, row_offset + i, col_offset*2, (void*)submatrix[i]);
    }

    return submatrix;
}

void recursive_complex_dmatmul(double **C, double **A, double **B, int M, int N, int P, int depth) {
    //printf("Recursive called at depth %d with M %d N %d P %d\n",depth,M,N,P);
    if (M <= TILE_SIZE && N <= TILE_SIZE && P <= TILE_SIZE) {
        base_complex_dmatmul(C, A, B, M, N, P);
        return;
    }

    const int MAX_DEPTH = 2;
    int m2 = M / 2;
    int n2 = N / 2;
    int p2 = P / 2;

    #pragma omp parallel if (depth == 0)
    {
        #pragma omp single
        {
            if (M >= N && M >= P) {
                //printf("CASE 1\n");
                double **C1 = create_complex_submatrix(C, 0, 0, m2, P);
                double **C2 = create_complex_submatrix(C, m2, 0, M - m2, P);
                double **A1 = create_complex_submatrix(A, 0, 0, m2, N);
                double **A2 = create_complex_submatrix(A, m2, 0, M - m2, N);

                #pragma omp task if (depth < MAX_DEPTH)
                recursive_complex_dmatmul(C1, A1, B, m2, N, P, depth + 1);

                #pragma omp task if (depth < MAX_DEPTH)
                recursive_complex_dmatmul(C2, A2, B, M - m2, N, P, depth + 1);

                #pragma omp taskwait

                free(C1);
                free(C2);
                free(A1);
                free(A2);
            } else if (P >= M && P >= N) {
                //printf("CASE 2\n");
                double **C1 = create_complex_submatrix(C, 0, 0, M, p2);
                double **C2 = create_complex_submatrix(C, 0, p2, M, P - p2);
                double **B1 = create_complex_submatrix(B, 0, 0, N, p2);
                double **B2 = create_complex_submatrix(B, 0, p2, N, P - p2);

                #pragma omp task if (depth < MAX_DEPTH)
                recursive_complex_dmatmul(C1, A, B1, M, N, p2, depth + 1);

                #pragma omp task if (depth < MAX_DEPTH)
                recursive_complex_dmatmul(C2, A, B2, M, N, P - p2, depth + 1);

                #pragma omp taskwait

                free(C1);
                free(C2);
                free(B1);
                free(B2);
            } else {
                //printf("CASE 3\n");
                double **A1 = create_complex_submatrix(A, 0, 0, M, n2);
                double **A2 = create_complex_submatrix(A, 0, n2, M, N - n2);
                double **B1 = create_complex_submatrix(B, 0, 0, n2, P);
                double **B2 = create_complex_submatrix(B, n2, 0, N - n2, P);

                #pragma omp task if (depth < MAX_DEPTH)
                recursive_complex_dmatmul(C, A1, B1, M, n2, P, depth + 1);

                #pragma omp task if (depth < MAX_DEPTH)
                recursive_complex_dmatmul(C, A2, B2, M, N - n2, P, depth + 1);

                #pragma omp taskwait

                free(A1);
                free(A2);
                free(B1);
                free(B2);
            }
        }
    }
}



// Function to create a row pointer array for complex matrix stored in contiguous memory
template <typename T>
std::complex<T> **create_row_pointers(std::complex<T> *matrix, int rows, int cols) {

    std::complex<T> **row_ptrs = (std::complex<T> **) malloc((rows + 1) * sizeof(std::complex<T> *));
    if (row_ptrs == NULL) {
        perror("Failed to allocate row pointers");
        exit(EXIT_FAILURE);
    }
    
    // Each row starts at 2 * (i * cols) + 1 (1-based indexing for complex) // TODO: fix
    for (int i = 1; i <= rows; i++) {
        row_ptrs[i] = &matrix[(i - 1) * cols];
    }
    return row_ptrs;
}



template <typename T>
void matrix_mult(std::complex<T> *C, std::complex<T> *A, std::complex<T> *B, int M, int N, int P)
{
    std::complex<T> **row_ptrs_C = create_row_pointers(C, M, P);
    std::complex<T> **row_ptrs_A = create_row_pointers(A, M, N);
    std::complex<T> **row_ptrs_B = create_row_pointers(B, N, P);

    recursive_matmul(row_ptrs_C, row_ptrs_A, row_ptrs_B, M, N, P, 0);

    free(row_ptrs_C);
    free(row_ptrs_A);
    free(row_ptrs_B);
};
