/*
 * MIT STAR Lab
 * Pure software implementations of Matrix Fourier Transform to help test hardware kernels (C++)
 * 
 * C code inherted from Nicholas Belsten.
 * Last modified on Feb 23, 2025 
 * by M.Subhi Abo Rdan (msubhi_a@mit.edu)
 */

#include "mft.h"

#define LINKER_EXPLICIT_INSTANTIATION_FLOAT
#ifdef LINKER_EXPLICIT_INSTANTIATION_FLOAT
  template void mft_forward(std::complex<float> **wavefront, float psf_pixelscale_lamD, int npsf, int npix);
  template void mft_reverse(std::complex<float> **wavefront, float psf_pixelscale_lamD, int npsf, int npix);
  template void complex_mat_mul(std::complex<float> *C, std::complex<float> *A, std::complex<float> *B, int M, int N, int P);
#endif // LINKER_EXPLICIT_INSTANTIATION_FLOAT

#define LINKER_EXPLICIT_INSTANTIATION_DOUBLE
#ifdef LINKER_EXPLICIT_INSTANTIATION_DOUBLE
 template void mft_forward(std::complex<double> **wavefront, double psf_pixelscale_lamD, int npsf, int npix);
 template void mft_reverse(std::complex<double> **wavefront, double psf_pixelscale_lamD, int npsf, int npix);
 template void complex_mat_mul(std::complex<double> *C, std::complex<double> *A, std::complex<double> *B, int M, int N, int P);
#endif // LINKER_EXPLICIT_INSTANTIATION_DOUBLE


template <typename T>
void mft_forward(std::complex<T> **wavefront, T psf_pixelscale_lamD, int npsf, int npix)
{

    T dx = 1.0 / (T)(npix);
    T du = psf_pixelscale_lamD;

    T* Xs = (T*) malloc(npix * sizeof(T));
    T* Us = (T*) malloc(npsf * sizeof(T));
    if (!Xs || !Us) {
        free(Xs);
        free(Us);
        return; // Handle allocation failure
    }

    for (int i = 0; i < npix; i++) {
        Xs[i] = ((T)(i) - (T)(npix) / 2.0) * dx;
    }

    for (int i = 0; i < npsf; i++) {
        Us[i] = ((T)(i) - (T)(npsf) / 2.0) * du;
    }

    std::complex<T>* Mx = (std::complex<T>*) malloc(npsf * npix * sizeof(std::complex<T>));
    std::complex<T>* My = (std::complex<T>*) malloc(npsf * npix * sizeof(std::complex<T>));
    if (!Mx || !My) {
        free(Mx);
        free(My);
        free(Xs);
        free(Us);
        return;
    }

    std::complex<T> I(0, 1);

    for (int i = 0; i < npsf; i++) {
        for (int j = 0; j < npix; j++) {
            int index = i * npix + j;
            T x_outer = Us[i] * Xs[j];
            Mx[index] = std::exp(-2 * static_cast<T>(M_PI) * I * x_outer) / (T)(npix);
        }
    }

    for (int i = 0; i < npix; i++) {
        for (int j = 0; j < npsf; j++) {
            int index = i * npsf + j;
            T y_outer = Xs[i] * Us[j];
            My[index] = std::exp(-2 * static_cast<T>(M_PI) * I * y_outer) * du;
        }
    }

    // Allocate temp storage
    std::complex<T>* temp = (std::complex<T>*) calloc(npsf * npix, sizeof(std::complex<T>));
    if (!temp) {
        free(Mx);
        free(My);
        free(Xs);
        free(Us);
        return;
    }

    // Perform matrix multiplications Mx@pupil@My
    complex_mat_mul(temp, *wavefront, My, npix, npix, npsf);

    // Free original wavefront before reallocation
    free(*wavefront);
    *wavefront = (std::complex<T>*) calloc(npsf * npsf, sizeof(std::complex<T>));

    complex_mat_mul(*wavefront, Mx, temp, npsf, npix, npsf);

    // Free memory
    free(temp);
    free(Mx);
    free(My);
    free(Xs);
    free(Us);
}


 // Note: The wavefront is deallocated then reallocated at a different size in this function
 template <typename T>
 void mft_reverse(std::complex<T> **wavefront, T psf_pixelscale_lamD, int npsf, int npix)
 {
    T dx = 1.0 / (T)(npix);
    T du = psf_pixelscale_lamD;

    T* Xs = (T*) malloc(npix * sizeof(T));
    T* Us = (T*) malloc(npsf * sizeof(T));
    if (!Xs || !Us) {
        free(Xs);
        free(Us);
        return; // Handle allocation failure
    }

    for (int i = 0; i < npix; i++) {
        Xs[i] = ((T)(i) - (T)(npix) / 2.0) * dx;
    }

    for (int i = 0; i < npsf; i++) {
        Us[i] = ((T)(i) - (T)(npsf) / 2.0) * du;
    }

    std::complex<T>* Mx = (std::complex<T>*) malloc(npsf * npix * sizeof(std::complex<T>));
    std::complex<T>* My = (std::complex<T>*) malloc(npsf * npix * sizeof(std::complex<T>));
    if (!Mx || !My) {
        free(Mx);
        free(My);
        free(Xs);
        free(Us);
        return;
    }

    std::complex<T> I(0, 1);


    for (int i = 0; i < npsf; i++) {
        for (int j = 0; j < npix; j++) {
            int index = i * npix + j;
            T yv = Us[i] * Xs[j];
            My[index] = std::exp(-2 * static_cast<T>(M_PI) * I * yv) / (T)npix;
        }
    }

    for (int i = 0; i < npix; i++) {
        for (int j = 0; j < npsf; j++) {
            int index = i * npsf + j;
            T ux = Xs[i] * Us[j];
            Mx[index] = std::exp(-2 * static_cast<T>(M_PI) * I * ux) * du;
        }
    }

    //return Mx@fpwf@My
    std::complex<T>* temp = (std::complex<T>*) calloc(npsf * npix, sizeof(std::complex<T>));
    if (!temp) {
        free(Mx);
        free(My);
        free(Xs);
        free(Us);
        return;
    }

    // Perform matrix multiplications Mx@pupil@My
    complex_mat_mul(temp, *wavefront, My, npsf, npsf, npix);
    complex_mat_mul(*wavefront,Mx,temp, npix, npsf, npix);

    free(temp);
    free(Mx);
    free(My);
    free(Xs);
    free(Us);
}









/********************  Matrix Complex Multiplication **********************/

#define TILE_SIZE 16  // Set tile size as a constant

template <typename T>
void print_matrix(std::complex<T> *mat, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << mat[i * cols + j] << "\t";
        }
        std::cout << std::endl;
    }
}

template <typename T>
void base_complex_dmatmul(std::complex<T> **C, std::complex<T> **A, std::complex<T> **B, int M, int N, int P)
{
    int m, n, p;
    for (m = 0; m < M; m++) {
        for (p = 0; p < P; p++) {
            std::complex<T> sum(0,0);

            for (n = 0; n < N; n++) {
                sum += A[m][n] * B[n][p];
            }
            // Accumulate on existing values in C
            C[m][p] += sum;
        }
    }
}

// Helper to create a submatrix pointer array for complex matrices
template <typename T>
std::complex<T> **create_complex_submatrix(std::complex<T> **matrix, int row_offset, int col_offset, int rows, int cols) 
{
    std::complex<T> **submatrix = (std::complex<T> **)malloc(rows * sizeof(std::complex<T> *));
    if (submatrix == NULL) {
        perror("Failed to allocate submatrix pointers");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; i++) {
        submatrix[i] = &matrix[row_offset + i][col_offset];
    }

    return submatrix;
}

template <typename T>
void recursive_complex_dmatmul(std::complex<T> **C, std::complex<T> **A, std::complex<T> **B, int M, int N, int P, int depth) 
{
    //printf("Recursive called at depth %d with M %d N %d P %d\n",depth,M,N,P);
    if (M <= TILE_SIZE && N <= TILE_SIZE && P <= TILE_SIZE) {
        base_complex_dmatmul(C, A, B, M, N, P);
        // std::cout<< "1.1 passed" << std::endl;
        return;
    }
    // std::cout<< "1 bug" << std::endl;

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
                std::complex<T> **C1 = create_complex_submatrix(C, 0, 0, m2, P);
                std::complex<T> **C2 = create_complex_submatrix(C, m2, 0, M - m2, P);
                std::complex<T> **A1 = create_complex_submatrix(A, 0, 0, m2, N);
                std::complex<T> **A2 = create_complex_submatrix(A, m2, 0, M - m2, N);

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
                std::complex<T> **C1 = create_complex_submatrix(C, 0, 0, M, p2);
                std::complex<T> **C2 = create_complex_submatrix(C, 0, p2, M, P - p2);
                std::complex<T> **B1 = create_complex_submatrix(B, 0, 0, N, p2);
                std::complex<T> **B2 = create_complex_submatrix(B, 0, p2, N, P - p2);

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
                std::complex<T> **A1 = create_complex_submatrix(A, 0, 0, M, n2);
                std::complex<T> **A2 = create_complex_submatrix(A, 0, n2, M, N - n2);
                std::complex<T> **B1 = create_complex_submatrix(B, 0, 0, n2, P);
                std::complex<T> **B2 = create_complex_submatrix(B, n2, 0, N - n2, P);

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
std::complex<T> **create_row_pointers(std::complex<T>  *matrix, int rows, int cols) 
{

    std::complex<T> **row_ptrs = (std::complex<T> **)malloc(rows * sizeof(std::complex<T>*));
    if (row_ptrs == NULL) {
        perror("Failed to allocate row pointers");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; i++) {
        row_ptrs[i] = &matrix[i*cols];
    }

    return row_ptrs;
}


template <typename T>
void complex_mat_mul(std::complex<T> *C, std::complex<T> *A, std::complex<T> *B, int M, int N, int P)
{
    std::complex<T> **row_ptrs_C = create_row_pointers(C, M, P);
    std::complex<T> **row_ptrs_A = create_row_pointers(A, M, N);
    std::complex<T> **row_ptrs_B = create_row_pointers(B, N, P);

    // std::cout<< "1 passed" << std::endl;
    recursive_complex_dmatmul(row_ptrs_C, row_ptrs_A, row_ptrs_B, M, N, P, 0);

    // std::cout<< "2 passed" << std::endl;
    free(row_ptrs_C);
    free(row_ptrs_A);
    free(row_ptrs_B);
}






////////////////// wrappers used for python testing

#ifdef PYTHON_TESTING_FLOATS
extern "C" {
    void mft_forward_float_wrapper(std::complex<float> *wavefront, std::complex<float> *wavefront_output, float psf_pixelscale_lamD, int npsf, int npix){
        std::complex<float>** wavefront_pp = (std::complex<float>**) malloc(sizeof(std::complex<float>*));
        *wavefront_pp = (std::complex<float>*) malloc(npix * npix * sizeof(std::complex<float>));
        memcpy(*wavefront_pp, wavefront, npix * npix * sizeof(std::complex<float>));

        mft_forward(wavefront_pp, psf_pixelscale_lamD, npsf, npix);

        memcpy(wavefront_output, *wavefront_pp, npsf * npsf * sizeof(std::complex<float>));
        free(wavefront_pp);
    }

    void mft_reverse_float_wrapper(std::complex<float> *wavefront, std::complex<float> *wavefront_output, float psf_pixelscale_lamD, int npsf, int npix){
        std::complex<float>** wavefront_pp = (std::complex<float>**) malloc(sizeof(std::complex<float>*));
        *wavefront_pp = (std::complex<float>*) malloc(npix * npix * sizeof(std::complex<float>));
        memcpy(*wavefront_pp, wavefront, npix * npix * sizeof(std::complex<float>));

        mft_reverse(wavefront_pp, psf_pixelscale_lamD, npsf, npix);

        memcpy(wavefront_output, *wavefront_pp, npsf * npsf * sizeof(std::complex<float>));
        free(wavefront_pp);
    }
    void complex_mat_mul_float_wrapper(std::complex<float> *C, std::complex<float> *A, std::complex<float> *B, int M, int N, int P){
        complex_mat_mul(C, A, B, M, N, P);
    }
}
#endif // PYTHON_TESTING_FLOATS


#ifdef PYTHON_TESTING_DOUBLES
extern "C" {
    void mft_forward_float_wrapper(std::complex<double> *wavefront, std::complex<double> *wavefront_output, double psf_pixelscale_lamD, int npsf, int npix){
        std::complex<double>** wavefront_pp = (std::complex<double>**) malloc(sizeof(std::complex<double>*));
        *wavefront_pp = (std::complex<double>*) malloc(npix * npix * sizeof(std::complex<double>));
        memcpy(*wavefront_pp, wavefront, npix * npix * sizeof(std::complex<double>));

        mft_forward(wavefront_pp, psf_pixelscale_lamD, npsf, npix);

        memcpy(wavefront_output, *wavefront_pp, npsf * npsf * sizeof(std::complex<double>));
        free(wavefront_pp);
    }

    void mft_reverse_double_wrapper(std::complex<double> *wavefront, std::complex<double> *wavefront_output, double psf_pixelscale_lamD, int npsf, int npix){
        std::complex<double>** wavefront_pp = (std::complex<double>**) malloc(sizeof(std::complex<double>*));
        *wavefront_pp = (std::complex<double>*) malloc(npix * npix * sizeof(std::complex<double>));
        memcpy(*wavefront_pp, wavefront, npix * npix * sizeof(std::complex<double>));

        mft_reverse(wavefront_pp, psf_pixelscale_lamD, npsf, npix);

        memcpy(wavefront_output, *wavefront_pp, npsf * npsf * sizeof(std::complex<double>));
        free(wavefront_pp);
    }
    void complex_mat_mul_double_wrapper(std::complex<double> *C, std::complex<double> *A, std::complex<double> *B, int M, int N, int P){
        complex_mat_mul(C, A, B, M, N, P);
    }
}
#endif // PYTHON_TESTING_DOUBLES
