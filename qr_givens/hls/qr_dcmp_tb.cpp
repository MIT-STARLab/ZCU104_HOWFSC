

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "qr_dcmp.h"

using namespace std;
const double TOLERANCE = 1e-5*N*N;



// Structure to represent a matrix
typedef struct {
    int rows, cols;
    data_t* data;
} Matrix;



// Allocate memory for a matrix
Matrix create_matrix(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (data_t*)malloc(rows * cols * sizeof(data_t));
    return mat;
}

// Free matrix memory
void free_matrix(Matrix* mat) {
    free(mat->data);
    mat->data = NULL;
}

// Access matrix element (row-major order)
data_t get_element(const Matrix* mat, int row, int col) {
    return mat->data[row * mat->cols + col];
}

// Set matrix element
void set_element(Matrix* mat, int row, int col, data_t value) {
    mat->data[row * mat->cols + col] = value;
}

// Print matrix
void print_matrix(const Matrix* mat) {
    for (int i = 0; i < mat->rows; ++i) {
        for (int j = 0; j < mat->cols; ++j) {
            printf("%.8f ", get_element(mat, i, j));
        }
        printf("\n");
    }
}

// Print matrix
void print_raw_matrix(data_t* mat) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            printf("%.8f ", mat[(i*N)+j]);
        }
        printf("\n");
    }
}

// Create an identity matrix
Matrix create_identity_matrix(int size) {
    Matrix mat = create_matrix(size, size);
    for (int i = 0; i < size; ++i) {
        set_element(&mat, i, i, 1.0);
    }
    return mat;
}

// Helper function to compare two matrices
bool compare_matrices(const Matrix *mat1, const Matrix *mat2, const char *matrix_name) {
    if (mat1->rows != mat2->rows || mat1->cols != mat2->cols) {
        std::cerr << "Matrix dimension mismatch in " << matrix_name << " comparison!" << std::endl;
        return false;
    }

    for (int i = 0; i < mat1->rows; ++i) {
        for (int j = 0; j < mat1->cols; ++j) {
            double diff = std::fabs(mat1->data[i * mat1->cols + j] - mat2->data[i * mat2->cols + j]);
            if (diff > TOLERANCE) {
                std::cerr << "Mismatch found in " << matrix_name << " at (" << i << ", " << j
                          << "): Software = " << mat1->data[i * mat1->cols + j]
                          << ", Hardware = " << mat2->data[i * mat2->cols + j]
                          << ", Difference = " << diff << std::endl;
                return false;
            }
        }
    }

    return true;
}

// Matrix multiplication: C = A * B
Matrix multiply_matrices(const Matrix* A, const Matrix* B) {
    if (A->cols != B->rows) {
        fprintf(stderr, "Matrix dimension mismatch for multiplication\n");
        exit(EXIT_FAILURE);
    }

    Matrix C = create_matrix(A->rows, B->cols);
    for (int i = 0; i < A->rows; ++i) {
        for (int j = 0; j < B->cols; ++j) {
            data_t sum = 0.0;
            for (int k = 0; k < A->cols; ++k) {
                sum += get_element(A, i, k) * get_element(B, k, j);
            }
            set_element(&C, i, j, sum);
        }
    }
    return C;
}

// Transpose of a matrix
Matrix transpose_matrix(const Matrix* mat) {
    Matrix T = create_matrix(mat->cols, mat->rows);
    for (int i = 0; i < mat->rows; ++i) {
        for (int j = 0; j < mat->cols; ++j) {
            set_element(&T, j, i, get_element(mat, i, j));
        }
    }
    return T;
}

// Givens rotation
void givens_rotation(data_t a, data_t b, data_t* cos, data_t* sin) {
    data_t hypot = sqrt(a * a + b * b);
    *cos = a / hypot;
    *sin = -b / hypot;
}

// QR decomposition using Givens rotations
void qr_givens(Matrix* A, Matrix* QT, Matrix* R) {
    // *R = create_matrix(A->rows, A->cols);
    // *QT = create_identity_matrix(A->rows);

    // Copy A into R
    for (int i = 0; i < A->rows; ++i) {
        for (int j = 0; j < A->cols; ++j) {
            set_element(R, i, j, get_element(A, i, j));
        }
    }

    for (int i = 0; i < A->cols - 1; ++i) {
        for (int j = i + 1; j < A->rows; ++j) {
            data_t cos, sin;
            givens_rotation(get_element(R, i, i), get_element(R, j, i), &cos, &sin);

            // Apply Givens rotation to rows i and j of R and QT
            for (int k = 0; k < A->cols; ++k) {
                data_t tempR = get_element(R, i, k) * cos + get_element(R, j, k) * (-sin);
                set_element(R, j, k, get_element(R, i, k) * sin + get_element(R, j, k) * cos);
                set_element(R, i, k, tempR);

                data_t tempQT = get_element(QT, i, k) * cos + get_element(QT, j, k) * (-sin);
                set_element(QT, j, k, get_element(QT, i, k) * sin + get_element(QT, j, k) * cos);
                set_element(QT, i, k, tempQT);
            }
        }
    }
}

void random_data_generator_array(data_t *data) {
    mt19937 gen(1703);
    uniform_real_distribution<> dis(100.0, 1000.0);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            data[i*N+j] = dis(gen);
        }
    }
}

void eye(data_t* mat){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if(i!=j){
                mat[i*N+j]=0.0;                
            }else{
                mat[i*N+j]=1.0;
            }
        }
    }
}


int main() {
    data_t* software_A = new data_t[MAT_SIZE];
    data_t* software_QT = new data_t[MAT_SIZE];
    data_t* software_R = new data_t[MAT_SIZE];

    data_t* hardware_A = new data_t[MAT_SIZE];
    data_t* hardware_QT = new data_t[MAT_SIZE];
    data_t* hardware_R = new data_t[MAT_SIZE];

    cout << "Generate Testing Data using Pure software implementaion" << endl;
    /////////////   Generate Testing Data using Pure software implementaion  //////////////


    random_data_generator_array(software_A);
    // for(int i=0; i<MAT_SIZE; i++){
    //     software_A[i]=i+1;
    // }
    eye(software_QT);

    for(int i=0; i<MAT_SIZE; i++){
        software_R[i]=software_A[i];

        hardware_A[i]=software_A[i];
        hardware_R[i]=software_R[i];
        hardware_QT[i]=software_QT[i];
   }
    cout << "Testing data generated" << endl;


    cout << "Call QRDCMP software function" << endl;
    /////////////////  Call QR DCMP function   //////////////
    Matrix A, QT, R;
    A.data = software_A;
    A.cols = N;
    A.rows = N;


    QT.data = software_QT;
    QT.cols = N;
    QT.rows = N;

    R.data=software_R;
    R.cols=N;
    R.rows=N;


    // Perform QR decomposition
    qr_givens(&A, &QT, &R);

    // Print results
    printf("QT:\n");
    print_matrix(&QT);

    printf("R:\n");
    print_matrix(&R);

    // Verify orthogonality
    Matrix QT_transpose = transpose_matrix(&QT);
    Matrix QQT = multiply_matrices(&QT, &QT_transpose);
    Matrix QTQ = multiply_matrices(&QT_transpose, &QT);

    printf("Q Q**T:\n");
    print_matrix(&QQT);

    printf("Q**T Q:\n");
    print_matrix(&QTQ);


    cout << "Running hardware kernel function" << endl;
    /////////////////  Run kernel   //////////////    
    krnl_qr_dcmp(hardware_QT, hardware_R);
    Matrix hardware_QT_mat, hardware_R_mat;
    
    hardware_QT_mat.data = hardware_QT;
    hardware_QT_mat.rows = N;
    hardware_QT_mat.cols = N;
    
    hardware_R_mat.data = hardware_R;
    hardware_R_mat.rows = N;
    hardware_R_mat.cols = N;

    // Print results
    printf("QT:\n");
    print_matrix(&hardware_QT_mat);

    printf("R:\n");
    print_matrix(&hardware_R_mat);


    cout << "Validating Results" << endl;
    /////////////////  Validate Results   //////////////    
    bool passed = true;


    // Compare QT matrices
    if (!compare_matrices(&QT, &hardware_QT_mat, "QT")) {
        passed = false;
    }

    // Compare R matrices
    if (!compare_matrices(&R, &hardware_R_mat, "R")) {
        passed = false;
    }

    if (passed) {
        cout << "All tests passed!" << endl;
    } else{
        cout <<"Tests failed "  << endl;
    }

    return 0;
}
