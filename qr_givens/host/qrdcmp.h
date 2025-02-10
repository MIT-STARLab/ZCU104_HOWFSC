
#ifndef QRDCMP_SW_H
#define QRDCMP_SW_H


using namespace std;


//Passing a command line argument now to set the below
//#define N 64
//#define MAT_SIZE (N * N)


typedef float data_t;


// Structure to represent a matrix
typedef struct {
    int rows, cols;
    data_t* data;
} Matrix;


Matrix create_matrix(int rows, int cols);
void free_matrix(Matrix* mat);
data_t get_element(const Matrix* mat, int row, int col);
void set_element(Matrix* mat, int row, int col, data_t value);
void print_matrix(const Matrix* mat);
void print_raw_matrix(data_t* mat, int n);
Matrix create_identity_matrix(int size);
void random_data_generator_array(data_t *data, int n);
void eye(data_t* mat, int n);
void qr_givens(Matrix* A, Matrix* QT, Matrix* R);

// Top Level QR DCMP Kernel
extern "C" void krnl_qr_dcmp(data_t *QT, data_t *R);



#endif //QRDCMP_SW_H