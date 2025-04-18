#include <iostream>
#include <cstdlib>
#include <cmath>
#include <stdexcept>

#include <cstring>     // Add for fabs
#include <malloc.h> 
#include <time.h>


#include "mxm_r.h" // Assume this includes the kernel declaration

const int M = MATRIX_M, N = MATRIX_N, P = MATRIX_P;
const int tile_size = TILE_SIZE;


void base_matmul(data_t* C, data_t* A, data_t* B,
                 int a_i_start, int k_start, int b_j_start,
                 int M, int N, int P) {
    // Calculate global tile boundaries based on current submatrix size
    const int i_end = a_i_start + TILE_SIZE;
    const int j_end = b_j_start + TILE_SIZE;
    const int k_end = k_start + TILE_SIZE;

    for (int i = a_i_start; i < i_end; i++) {
        for (int j = b_j_start; j < j_end; j++) {
            data_t acc = 0.0f;
            for (int k = k_start; k < k_end; k++) {
                acc += A[i * N + k] * B[k * P + j];
            }
            C[i * P + j] += acc;
        }
    }
}

void recursive_matmul(data_t* C, data_t* A, data_t* B,
                     int a_i_start, int k_start, int b_j_start,
                     int curr_M, int curr_N, int curr_P,
                     int M, int N, int P, int depth) {
    
    if (curr_M <= TILE_SIZE && curr_N <= TILE_SIZE && curr_P <= TILE_SIZE) {
        base_matmul(C, A, B, a_i_start, k_start, b_j_start, M, N, P);
        return;
    }
    {
        // Modified splitting logic to prioritize M/P over N
        if (curr_M > TILE_SIZE) {
            int split = curr_M / 2;
            recursive_matmul(C, A, B, a_i_start, k_start, b_j_start, split, curr_N, curr_P, M, N, P, depth + 1);
            recursive_matmul(C, A, B, a_i_start + split, k_start, b_j_start, curr_M - split, curr_N, curr_P, M, N, P, depth + 1);
        }
        else if (curr_P > TILE_SIZE) {
            int split = curr_P / 2;
            recursive_matmul(C, A, B, a_i_start, k_start, b_j_start, curr_M, curr_N, split, M, N, P, depth + 1);
            recursive_matmul(C, A, B, a_i_start, k_start, b_j_start + split, curr_M, curr_N, curr_P - split, M, N, P, depth + 1);
        }
        else { // Only split N if M/P are already at tile size
            int split = curr_N / 2;
            recursive_matmul(C, A, B, a_i_start, k_start, b_j_start, curr_M, split, curr_P, M, N, P, depth + 1);
            recursive_matmul(C, A, B, a_i_start, k_start + split, b_j_start, curr_M, curr_N - split, curr_P, M, N, P, depth + 1);
        }
    }
}

void recursive_matmul_hls(data_t* C, data_t* A, data_t* B,
                     int a_i_start, int k_start, int b_j_start,
                     int curr_M, int curr_N, int curr_P,
                     int M, int N, int P, int depth) {
    
    if (curr_M <= TILE_SIZE && curr_N <= TILE_SIZE && curr_P <= TILE_SIZE) {
        krnl_mxm_r(C, A, B, a_i_start, k_start, b_j_start, M, N, P);
        return;
    }
    {
        // Modified splitting logic to prioritize M/P over N
        if (curr_M > TILE_SIZE) {
            int split = curr_M / 2;
            recursive_matmul_hls(C, A, B, a_i_start, k_start, b_j_start, split, curr_N, curr_P, M, N, P, depth + 1);
            recursive_matmul_hls(C, A, B, a_i_start + split, k_start, b_j_start, curr_M - split, curr_N, curr_P, M, N, P, depth + 1);
        }
        else if (curr_P > TILE_SIZE) {
            int split = curr_P / 2;
            recursive_matmul_hls(C, A, B, a_i_start, k_start, b_j_start, curr_M, curr_N, split, M, N, P, depth + 1);
            recursive_matmul_hls(C, A, B, a_i_start, k_start, b_j_start + split, curr_M, curr_N, curr_P - split, M, N, P, depth + 1);
        }
        else { // Only split N if M/P are already at tile size
            int split = curr_N / 2;
            recursive_matmul_hls(C, A, B, a_i_start, k_start, b_j_start, curr_M, split, curr_P, M, N, P, depth + 1);
            recursive_matmul_hls(C, A, B, a_i_start, k_start + split, b_j_start, curr_M, curr_N - split, curr_P, M, N, P, depth + 1);
        }
    }
}



data_t* allocate_matrix(int rows, int cols) {
    // Correct parameter order: size first, alignment second
    data_t* mat = (data_t*) aligned_alloc(rows * cols * sizeof(data_t), 64);
    if (!mat) {
        std::cerr << "Memory allocation failed!" << std::endl;
        exit(1);
    }
    return mat;
}

// In free_matrix function
void free_matrix(data_t* matrix) {
    free(matrix);  // Use aligned free for aligned allocations
}

void set_matrix(data_t* matrix, int rows, int cols) {
    memset(matrix, 0, rows * cols * sizeof(data_t));
}


// ... [keep configure_omp, allocate_matrix, free_matrix, zero_matrix unchanged]
void standard_matmul(data_t* C, data_t* A, data_t* B, int M, int N, int P) {
    //set_matrix(C, M, P);
    recursive_matmul(C, A, B, 0, 0, 0, M, N, P, M, N, P, 0);
}

void standard_matmul_nr(data_t* C, data_t* A, data_t* B, int M, int N, int P) {
    // Loop over rows of A (and C)
    for (int i = 0; i < M; i++) {
        // Loop over columns of B (and C)
        for (int j = 0; j < P; j++) {
            data_t sum = 0;
            // Compute dot product of the i-th row of A and the j-th column of B
            for (int k = 0; k < N; k++) {
                sum += A[i * N + k] * B[k * P + j];
            }
            C[i * P + j] = sum;
        }
    }
}

// ... [keep configure_omp, allocate_matrix, free_matrix, zero_matrix unchanged]
void hls_matmul(data_t* C, data_t* A, data_t* B, int M, int N, int P) {
    //set_matrix(C, M, P);
    recursive_matmul_hls(C, A, B, 0, 0, 0, M, N, P, M, N, P, 0);
}


// Compare sub-block of C_hls with C_gold
bool compare_result( data_t* C_hls,  data_t* C_gold,
                      int M, int P,   // full dimensions of C matrices
                      float tolerance) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < P; j++) {
            data_t diff = std::fabs(C_hls[(i) * P + (j)]
                                    - C_gold[(i) * P + (j)]);
            if (diff > tolerance) {
                std::cerr << "Mismatch at ("
                          << i << ", " << j << "): "
                          << "C_hls = " << C_hls[(i) * P + (j)]
                          << ", C_gold = " << C_gold[(i) * P + (j)]
                          << std::endl;
                return false;
            }
        }
    }
    return true;
}

void print_matrix(data_t* mat, int n, int cols){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << mat[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }

}

void initialize_matrix_rand(data_t *matrix, int rows, int cols, int minVal, int maxVal) {
    // Seed the random number generator (typically done once)
    srand(time(NULL));
    
    //int totalElements = rows * cols;
    int range = maxVal - minVal + 1;
    
    for (int i = 0; i < rows*cols; i++) {
        int random = (rand() % range) + minVal;
        matrix[i] = random;
    }
}

void initialize_matrix_one(data_t *matrix, int rows, int cols) {
    // Seed the random number generator (typically done once)
    srand(time(NULL));
    
    //int totalElements = rows * cols;
    //int range = maxVal - minVal + 1;
    
    for (int i = 0; i < rows*cols; i++) {
        //int random = (rand() % range) + minVal;
        matrix[i] = 1.0f;
    }
}

int main() {
    try {
        // Full matrix dimensions (larger than the tile)


        // Allocate matrices as contiguous memory
        data_t* A = new data_t[M*N];
        data_t* B = new data_t[N*P];
        data_t* C_gold = new data_t[M*P];
        data_t* C_hls = new data_t[M*P];

    
        // Initialize A and B
        std::cout << "Initializing matrices..." << std::endl;

        //for (int i = 0; i < M*N; i++) A[i] = 1.0f;
        //for (int i = 0; i < N*P; i++) B[i] = 1.0f;
        initialize_matrix_one(A, M, N);
        initialize_matrix_one(B, N, P);

        for (int i = 0; i < M*P; i++) C_hls[i] = 0.0f;
        for (int i = 0; i < M*P; i++) C_gold[i] = 0.0f;


        std::cout << "A:" << std::endl;
        print_matrix(A, 16, MATRIX_P);
        std::cout << "B:" << std::endl;
        print_matrix(B, 16, MATRIX_P);
        std::cout << "C_gold:" << std::endl;
        print_matrix(C_gold, 16, MATRIX_P);
        std::cout << "C_hls:" << std::endl;
        print_matrix(C_hls, 16, MATRIX_P);

        // Compute the golden result for the entire large matrix
        std::cout << "Computing golden result..." << std::endl;
        standard_matmul_nr(C_gold, A, B, M, N, P);

        // Run the kernel on the n x n sub-block
        std::cout << "Running HLS kernel on sub-block..." << std::endl;
        // krnl_mxm_r(
        //     A, B, C_hls, 
        //     start_i, start_j,
        //     M, N, P
        // );
        hls_matmul(C_hls, A, B, M, N, P);

        

        // Now print the first 4x4 sub-block
        std::cout << "Printing first 4x4 sub-block of C_gold and C_hls for inspection:\n";
        std::cout << "C_gold (first 4x4):" << std::endl;
        print_matrix(C_gold, 16, MATRIX_P);

        std::cout << "C_hls (first 4x4):" << std::endl;
        print_matrix(C_hls, 16, MATRIX_P);

        // Compare only the n x n sub-block
        std::cout << "Comparing tile results..." << std::endl;
        float tolerance = 1e-3;
        bool is_equal = compare_result(C_hls, C_gold, M, P, tolerance);

        if (is_equal) {
            std::cout << "Test PASSED: Sub-block matches golden result." << std::endl;
        } else {
            std::cerr << "Test FAILED: Sub-block does not match golden result." << std::endl;
        }

        // Free allocated memory
        // free_matrix(A);
        // free_matrix(B);
        // free_matrix(C_gold);
        // free_matrix(C_hls);

        delete[] A;
        delete[] B;
        delete[] C_gold;
        delete[] C_hls;


    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
