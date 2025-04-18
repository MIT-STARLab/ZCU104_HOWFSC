#include <cstdio>
#include <sys/syscall.h>      /* Definition of SYS_* constants */
#include <time.h>

#include <stdlib.h>
#include <iostream>
#include <string>

#include <complex>
#include <cmath>
#include <algorithm>

// XRT includes
#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

//#define SEQUENTIAL_TIME_TEST
//#define RECURSIVE_TIME_TEST

using namespace std;

typedef float data_t;

//These only affect the test bench operations
// #define MATRIX_M 1024
// #define MATRIX_N 1024
// #define MATRIX_P 1024

#define TILE_SIZE 128
// #define MAT_A_SIZE MATRIX_M*MATRIX_N
// #define MAT_B_SIZE MATRIX_N*MATRIX_P
// #define MAT_C_SIZE MATRIX_M*MATRIX_P

/////// Timing Helpers ///////
#define TIMER(label)  timespec label; syscall(SYS_clock_gettime, CLOCK_MONOTONIC, &label)
#define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0

/////// Error Profiling Helpers ///////
static const char*    STR_PASSED  = "PASSED:    ";
static const char*    STR_FAILED   = "FAILED: ";
static const char*    STR_INFO    = "INFO:  ";
static const char*    STR_USAGE   = "USAGE: ";
static const char*    STR_RESULTS   = "RESULTS: ";

static const char* error_message =
    "Error:\tResults Mismatch:\n"
    "\tIndex i = %d; CPU Results (%f, %f); Device Result (%f, %f)\n"
    "\tRelative error in real component = %f\n"
    "\tRelative error in imag component = %f\n";

static const char* timing_summary = 
    "Matrix Multiplication Run Time:\n"
    "-------------------------\n"
    "Matrix Dimension Size (M,N.P) =   (%d, %d, %d)\n"
    "Total CPU Time        =   %f\n"
    "Total FPGA Time       =   %f\n"
    "->Input Load Time     =   %f\n"
    "->Excution Time       =   %f\n"
    "->Output Load Time    =   %f\n";

static const char* final_results =
    "\n**********Timing Summary*********\n"
    "Kernel name                        = Matrix Multiplication\n"
    "Input Data Size                    = (%d, %d)\n"
    "Number of Run Rounds               = %d\n"
    "FPGA Average Run Time (With Sync)  = %f\n"
    "FPGA Average Run Time (No   Sync)  = %f\n"
    "CPU  Average Run Time              = %f\n"
    "**********************\n";



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

void recursive_matmul_hls(xrt::kernel &krnl,
                        xrt::bo &bo_C, xrt::bo &bo_A, xrt::bo &bo_B,
                        int a_i_start, int k_start, int b_j_start,
                        int curr_M, int curr_N, int curr_P,
                        int M, int N, int P, int depth) {
    
    if (curr_M <= TILE_SIZE && curr_N <= TILE_SIZE && curr_P <= TILE_SIZE) {
        //TIMER(startKernel);
        auto run = krnl(bo_C, bo_A, bo_B, a_i_start, k_start, b_j_start, M, N, P); //We need i, j, and k here. Must specify sub blocks in all three matrices
        run.wait();
        //TIMER(endKernel);
        //std::cout << STR_INFO << "TIME KERNEL " << ELAPSED(endKernel, startKernel) << std::endl;
        return;
    }
    {
        // Modified splitting logic to prioritize M/P over N
        if (curr_M > TILE_SIZE) {
            int split = curr_M / 2;
            recursive_matmul_hls(krnl, bo_C, bo_A, bo_B, a_i_start, k_start, b_j_start, split, curr_N, curr_P, M, N, P, depth + 1);
            recursive_matmul_hls(krnl, bo_C, bo_A, bo_B, a_i_start + split, k_start, b_j_start, curr_M - split, curr_N, curr_P, M, N, P, depth + 1);
        }
        else if (curr_P > TILE_SIZE) {
            int split = curr_P / 2;
            recursive_matmul_hls(krnl, bo_C, bo_A, bo_B, a_i_start, k_start, b_j_start, curr_M, curr_N, split, M, N, P, depth + 1);
            recursive_matmul_hls(krnl, bo_C, bo_A, bo_B, a_i_start, k_start, b_j_start + split, curr_M, curr_N, curr_P - split, M, N, P, depth + 1);
        }
        else { // Only split N if M/P are already at tile size
            int split = curr_N / 2;
            recursive_matmul_hls(krnl, bo_C, bo_A, bo_B, a_i_start, k_start, b_j_start, curr_M, split, curr_P, M, N, P, depth + 1);
            recursive_matmul_hls(krnl, bo_C, bo_A, bo_B, a_i_start, k_start + split, b_j_start, curr_M, curr_N - split, curr_P, M, N, P, depth + 1);
        }
    }
}

// void recursive_matmul_hls(
//     xrt::kernel &krnl,
//     xrt::bo &boC, xrt::bo &boA, xrt::bo &boB,
//     int i_start, int j_start,
//     int M, int P,
//     int full_M, int full_N, int full_P)
// {
//     // Base case if M and P fit in tile
//     if (M <= TILE_SIZE && P <= TILE_SIZE)
//     {
//         // Now loop over N in increments (or all at once if your kernel can handle it).
//         // E.g., do partial sums in the kernel for sub-ranges of k.
//         for (int k = 0; k < full_N; k += TILE_SIZE) {
//             int k_block_size = std::min(TILE_SIZE, full_N - k);

//             // Launch kernel to do partial sum:
//             //   C(i_start..i_start+M, j_start..j_start+P) +=
//             //        A(i_start.., k_block_start..) × B(k_block_start.., j_start..);
//             auto run = krnl(boC, boA, boB,
//                             i_start, k, j_start, 
//                             M, N, P);
//             run.wait(); // block until done
//         }
//         return;
//     }

//     // Recursive splitting of M or P
//     if (M > TILE_SIZE) {
//         int half = M / 2;
//         recursive_matmul_hls(krnl, boC, boA, boB, i_start, j_start,         half, P, full_M, full_N, full_P);
//         recursive_matmul_hls(krnl, boC, boA, boB, i_start + half, j_start, M - half, P, full_M, full_N, full_P);
//     }
//     else if (P > TILE_SIZE) {
//         int half = P / 2;
//         recursive_matmul_hls(krnl, boC, boA, boB, i_start,         j_start, M, half,     full_M, full_N, full_P);
//         recursive_matmul_hls(krnl, boC, boA, boB, i_start, j_start + half, M, P - half, full_M, full_N, full_P);
//     }
// }



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

// Initialize a matrix with a fixed value (in this example: 1.0f).
void initialize_matrix(data_t* matrix, int rows, int cols, data_t value) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i * cols + j] = value; 
        }
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



void run_with_recursion(
    xrt::kernel &krnl,
    xrt::bo &bo_C, xrt::bo &bo_A, xrt::bo &bo_B,
    int M, int N, int P) {
    
    // Synchronize buffer content with device side
    bo_A.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_B.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    // Initiate the recursion
   // TIMER(startrec);
    recursive_matmul_hls(krnl, bo_C, bo_A, bo_B, 0, 0, 0, M, N, P, M, N, P, 0);
    //TIMER(endrec);
    //std::cout << STR_INFO << "TIME WITHOUT SYNC " << ELAPSED(endrec, startrec) << std::endl;

    // Retrieving the results
    bo_C.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
}

int main(int argc, char** argv) {

    // Check command line arguments
    if(argc < 2) {
        std::cout << STR_USAGE << argv[0] <<" <binary_container.xclbin> <rounds>" <<std::endl;
        return EXIT_FAILURE;
    }

    ////////////////     Initialize XRT device and load xclbin    //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    // Loading the xclbin
    char* xclbinFilename = argv[1];
    int rounds = (argc < 3)? 10: std::stoi(argv[2]);
    int multiple_m = (argc < 4)? 1: std::stoi(argv[3]);
    int multiple_n = (argc < 5)? 1: std::stoi(argv[4]);
    int multiple_p = (argc < 6)? 1: std::stoi(argv[5]);

    //open the device and return its handle
    unsigned int dev_index = 0;
    auto my_device = xrt::device(dev_index);
    std::cout << STR_PASSED << "auto my_device = xrt::device(" << dev_index << ")" << std::endl;


    //load the xclbin file from disk. This will return the UUID of the xclbin.
    auto xclbin_uuid = my_device.load_xclbin(xclbinFilename);
    std::cout << STR_PASSED << "auto xclbin_uuid = my_device.load_xclbin(" << xclbinFilename << ")" << std::endl;



    ////////////////        Set Kernel and argument               //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    int MATRIX_M = TILE_SIZE * multiple_m;
    int MATRIX_N = TILE_SIZE * multiple_n;
    int MATRIX_P = TILE_SIZE * multiple_p;

    int MAT_A_SIZE = MATRIX_M * MATRIX_N;
    int MAT_B_SIZE = MATRIX_N * MATRIX_P;
    int MAT_C_SIZE = MATRIX_M * MATRIX_P;

    // Starting the (HLS) PL kernel
    auto krnl = xrt::kernel(my_device, xclbin_uuid, "krnl_mxm_r");

    std::cout << STR_PASSED << "auto krnl = xrt::kernel(my_device, xclbin_uuid, \"krnl_mxm_r:{krnl_mxm_r_1}\")" << std::endl;

    std::cout << STR_INFO << "Allocate Buffers in Global Memory" << std::endl;

    // Compute the size of arrays in bytes
    size_t size_a_in_bytes = MAT_A_SIZE * sizeof(float);
    size_t size_b_in_bytes = MAT_B_SIZE * sizeof(float);
    size_t size_c_in_bytes = MAT_C_SIZE * sizeof(float);

    std::cout << STR_INFO << "Matrix A size in words  = " << MAT_A_SIZE  * 2   << std::endl;
    std::cout << STR_INFO << "Matrix A size in bytes  = " << size_a_in_bytes    << std::endl;

    std::cout << STR_INFO << "Matrix A size in words  = " << MAT_B_SIZE  * 2   << std::endl;
    std::cout << STR_INFO << "Matrix A size in bytes  = " << size_b_in_bytes    << std::endl;

    std::cout << STR_INFO << "Matrix A size in words  = " << MAT_C_SIZE  * 2   << std::endl;
    std::cout << STR_INFO << "Matrix A size in bytes  = " << size_c_in_bytes    << std::endl;
                                                                                              
    auto bo_A  = xrt::bo(my_device, size_a_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1));      //kernel argument 0: A matrix array
    auto bo_B  = xrt::bo(my_device, size_b_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2));      //kernel argument 1: B matrix array
    auto bo_C = xrt::bo(my_device, size_c_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(0));      //kernel argument 2: C/output matrix array


    std::cout << STR_PASSED << "auto bo_A  = xrt::bo(my_device, size_a_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1)) (=" << krnl.group_id(0) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_B   = xrt::bo(my_device, size_b_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2)) (=" << krnl.group_id(1) << "))" << std::endl;
    std::cout << STR_PASSED << "auto bo_C = xrt::bo(my_device, size_c_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(0)) (=" << krnl.group_id(2) << "))" << std::endl;

    //Map the contents of the buffer object into host memory
    auto bo_A_map  = bo_A.map<data_t*>();
    auto bo_B_map  = bo_B.map<data_t*>();
    auto bo_C_map  = bo_C.map<data_t*>();

    std::cout << STR_PASSED << "auto bo_A_map  = bo_A.map<data_t>()" << std::endl;
    std::cout << STR_PASSED << "auto bo_B_map  = bo_B.map<data_t>()" << std::endl;
    std::cout << STR_PASSED << "auto bo_C_map  = bo_C.map<data_t>()" << std::endl;




    ////////////////              Generate First Testing Data            /////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << STR_INFO << "Generating Testing Data..." << std::endl;

    data_t* A_software = new data_t[MAT_A_SIZE];
    data_t* B_software = new data_t[MAT_B_SIZE];
    data_t* C_software = new data_t[MAT_C_SIZE];
    

    initialize_matrix(A_software, MATRIX_M, MATRIX_N, (data_t) 1.0f);
    initialize_matrix(B_software, MATRIX_N, MATRIX_P, (data_t) 1.0f);
    initialize_matrix(C_software, MATRIX_M, MATRIX_P, (data_t) 0.0f);

    std::cout << STR_PASSED << "Testing Data Generated" << std::endl;
    
    std::cout << STR_INFO << "Filling Argument Buffers with input data" << std::endl;
    data_t fill_val = (data_t)0.00f;
    fill(bo_A_map, bo_A_map + MAT_A_SIZE, fill_val);
    fill(bo_B_map, bo_B_map + MAT_B_SIZE, fill_val);
    fill(bo_C_map, bo_C_map + MAT_C_SIZE, fill_val);

    // Fill the data in the buffers
    std::cout << "Filling Buffers with input data" << std::endl;
    for (int i = 0; i < MAT_A_SIZE; i++) {
        bo_A_map[i] = A_software[i];
    }
    for (int i = 0; i < MAT_B_SIZE; i++) {
        bo_B_map[i] = B_software[i];
    }

  
    ////////////////              CPU Excution                 //////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////


    std::cout << STR_INFO << "Running on CPU First" << std::endl;
    TIMER(cpuTimeStart);
    standard_matmul_nr(C_software, A_software,  B_software, MATRIX_M, MATRIX_N, MATRIX_P);
    TIMER(cpuTimeEnd);
    std::cout << STR_PASSED << "CPU Done" << std::endl;
    


    ///////////////     Recursive Kernel Runs ala Recursive Call for timing      ////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    printf("Going to run with recursion\n");
    TIMER(recursiveStartTime);
    run_with_recursion(krnl, bo_C, bo_A, bo_B, MATRIX_M, MATRIX_N, MATRIX_P);
    //recursive_matmul_hls(krnl, bo_C, bo_A, bo_B, 0, 0, 0, MATRIX_M, MATRIX_N, MATRIX_P, MATRIX_M, MATRIX_N, MATRIX_P, 0);
    TIMER(recursiveEndTime);
    printf("Ran recursive of size %d in %f s\n", MATRIX_M, ELAPSED(recursiveEndTime, recursiveStartTime));

    
    //Comparing Results
    std::cout << "Comparing tile results..." << std::endl;
    float tolerance = 1e-3;
    bool is_equal = compare_result(bo_C_map, C_software, MATRIX_M, MATRIX_P, tolerance);

    std::cout << "Printing first 16x16 sub-blocks of C inspection:\n";
    std::cout << "C_cpu (first 4x4):" << std::endl;
    print_matrix(C_software, 4, MATRIX_P);

    std::cout << "C_hls (first 4x4):" << std::endl;
    print_matrix(bo_C_map, 4, MATRIX_P);

    if (is_equal) {
        std::cout << "Test PASSED: Sub-block matches golden result." << std::endl;
    } else {
        std::cerr << "Test FAILED: Sub-block does not match golden result." << std::endl;
    }




    ///////////////     Multiple Kernel Runs on Random Data for timing      ////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    std::cout << STR_INFO << "RUN KERNEL NUMBER OF TIMES =  " << rounds << std::endl;

    double device_exc_time_average = 0;
    double device_tot_time_average = 0;
    double cpu_time_average = 0;

    for (int r=0; r<rounds; r++){
        std::cout << STR_INFO << "STARTING ROUND : " << r << std::endl;


        TIMER(cpuTimeStart);
        //fft2d(software_generated_input_data, MAT_ROWS, MAT_COLS, invert, scale);
        //auto run = krnl(bo_A, bo_B, bo_C, start_i, start_j, M, N, P);

        TIMER(cpuTimeEnd);

        cpu_time_average += ELAPSED(cpuTimeEnd, cpuTimeStart)/rounds;


        TIMER(syncInputTime);
        bo_A.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        bo_B.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        fill(bo_C_map, bo_C_map + MAT_C_SIZE, fill_val);
        TIMER(syncInputTimeEnd);

        TIMER(runTime);
        recursive_matmul_hls(krnl, bo_C, bo_A, bo_B, 0, 0, 0, MATRIX_M, MATRIX_N, MATRIX_P, MATRIX_M, MATRIX_N, MATRIX_P, 0);
        TIMER(runTimeEnd);

        TIMER(syncOutputTime);
        bo_C.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
        TIMER(syncOutputEnd);
        
        //Comparing Results
        std::cout << "Comparing tile results..." << std::endl;
        float tolerance = 1e-3;
        bool is_equal = compare_result(bo_C_map, C_software, MATRIX_M, MATRIX_P, tolerance);

        std::cout << "Printing first 16x16 sub-blocks of C inspection:\n";
        std::cout << "C_cpu (first 4x4):" << std::endl;
        print_matrix(C_software, 4, MATRIX_P);

        std::cout << "C_hls (first 4x4):" << std::endl;
        print_matrix(bo_C_map, 4, MATRIX_P);

        if (is_equal) {
            std::cout << "Test PASSED: Sub-block matches golden result." << std::endl;
        } else {
            std::cerr << "Test FAILED: Sub-block does not match golden result." << std::endl;
        }


       // TIMING
        printf(timing_summary, MATRIX_M, MATRIX_N, MATRIX_P, 
                           ELAPSED(cpuTimeEnd, cpuTimeStart),
                           ELAPSED(syncOutputEnd, syncInputTime), 
                           ELAPSED(syncInputTimeEnd, syncInputTime), 
                           ELAPSED(runTimeEnd, runTime), 
                           ELAPSED(syncOutputEnd, syncOutputTime));
        
        device_exc_time_average += ELAPSED(runTimeEnd, runTime)/rounds;
        device_tot_time_average += ELAPSED(syncOutputEnd, syncInputTime)/rounds;
    }

    printf(final_results, MATRIX_M, MATRIX_P,
                          rounds,
                          device_tot_time_average,
                          device_exc_time_average,
                          cpu_time_average);



    delete[] A_software;
    delete[] B_software;
    delete[] C_software;

    return 0;


}