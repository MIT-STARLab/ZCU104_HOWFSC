#include "matrixvector.h"
#include <iostream>
#include <cstdlib> // For rand() and srand()
#include <ctime> // For time() to seed rand()
#include <cmath> // for fabs()
#include <iomanip>


// Simple function to compute matrix-vector multiplication
void simple_matrix_vector_mul(const data_t* A, const data_t* B, data_t* C, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        data_t tmp = 0.0;
        data_t elem = 0.0;
        C[i] =0; // Initialize the result element
        for (int j = 0; j < cols; j++) {
            data_t myA = A[i * cols + j];
            data_t myB = B[j];
            elem = A[i * cols + j] * B[j];
            C[i] += elem;
            tmp = C[i]; //just for debug
        }
    }
}

// Function to print a matrix stored in a 1D array
void print_matrix(const char* name, const data_t* A, int rows, int cols) {
    std::cout << "Matrix " << name << ":" << std::endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << std::setw(10) << A[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }
}

// Function to compare the results of the HLS kernel and the simple implementation
bool compare_results(const data_t* C1, const data_t* C2, int size, data_t tolerance = 1e-5) {
    for (int i = 0; i < size; ++i) {
        if (fabs(C1[i] - C2[i]) > tolerance) {
            return false; // Results differ
        }
    }
    return true; // Results are the same
}

int main() {
    srand(1); // Seed random number generator

    // Allocate memory for matrix A, vector B, and result vector C
    data_t* A = new data_t[MATRIX_ROWS * MATRIX_COLS];
    data_t* B = new data_t[MATRIX_COLS];
    data_t* C = new data_t[MATRIX_ROWS];

    // Initialize matrix A and vector B with random data
    for(int i = 0; i < MATRIX_ROWS * MATRIX_COLS; i++) {
        //A[i] = rand() % 10; // Random values between 0 and 9
        A[i] = i%10;
    }
    for(int i = 0; i < MATRIX_COLS; i++) {
        //B[i] = rand() % 10; // Random values between 0 and 9
        B[i] = i%10;

    }

    //print_matrix("A", A, MATRIX_ROWS, MATRIX_COLS);
    //print_matrix("B", B, MATRIX_COLS, 1);

    // Call the kernel function
    krnl_inner_product(C, A, B, MATRIX_ROWS, MATRIX_COLS);

    // Output the results
    std::cout << "Result vector C:" << std::endl;
    for(int i = 0; i < MATRIX_ROWS; i++) {
        std::cout << C[i] << std::endl;
    }

    
    // Allocate memory for the result of the simple function
    data_t* C_simple = new data_t[MATRIX_ROWS];

    // Perform the simple matrix-vector multiplication
    simple_matrix_vector_mul(A, B, C_simple, MATRIX_ROWS, MATRIX_COLS);

    // Output the simple results
    std::cout << "Simple result vector C:" << std::endl;
    for(int i = 0; i < MATRIX_ROWS; i++) {
        std::cout << C_simple[i] << std::endl;
    }

    // Compare the results
    if (compare_results(C, C_simple, MATRIX_ROWS)) {
        std::cout << "The results match." << std::endl;
    } else {
        std::cout << "The results do not match!" << std::endl;
    }

    // Clean up
    delete[] A;
    delete[] B;
    delete[] C;

    delete[] C_simple;

    return 0;
}
