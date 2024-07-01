#/*
#Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
#SPDX-License-Identifier: X11
#*/

//#include "cmdlineparser.h"
#include <iostream>
#include <cstring>
#include <thread>
#include <chrono>

// XRT includes
#include "xrt/xrt_bo.h"
#include <experimental/xrt_xclbin.h>
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"
// #include "experimental/xrt_queue.h"

using namespace std;

#define DATA_SIZE 7


int main(int argc, char** argv) {

    cout << "argc = " << argc << endl;
	for(int i=0; i < argc; i++){
	    cout << "argv[" << i << "] = " << argv[i] << endl;
	}

    // Read settings
    string binaryFile = argv[1];
    int device_index = 0;

    cout << "Open the device" << device_index << endl;
    auto device = xrt::device(device_index);
    cout << "Load the xclbin " << binaryFile << endl;
    auto uuid = device.load_xclbin(binaryFile);

    size_t vector_size_bytes = sizeof(int) * DATA_SIZE;

    auto krnl = xrt::kernel(device, uuid, "vadd_new", xrt::kernel::cu_access_mode::exclusive);

    std::cout << "Allocate Buffer in Global Memory\n";

    auto boIn1 = xrt::bo(device, vector_size_bytes, krnl.group_id(0));
    auto boIn2 = xrt::bo(device, vector_size_bytes, krnl.group_id(1));
    auto boOut = xrt::bo(device, vector_size_bytes, krnl.group_id(2));

    int input_buff1[DATA_SIZE] = {1,2,3,4,5,6,7};
    int input_buff2[DATA_SIZE] = {1,2,3,4,5,6,7};
    int bufReference[DATA_SIZE] = {2,4,6,8,10,12,14};
    int output_buff[DATA_SIZE];

    std::cout << "synchronize input buffer data to device global memory\n";
    boIn1.write(input_buff1);
    boIn1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    boIn2.write(input_buff2);
    boIn2.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    std::cout << "Execution of the kernel\n";
    auto run = krnl(boIn1, boIn2, boOut);

    run.wait();

    boOut.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    boOut.read(output_buff);


    // // Map the contents of the buffer object into host memory
    // auto bo0_map = boIn1.map<int*>();
    // auto bo1_map = boIn2.map<int*>();
    // auto bo2_map = boOut.map<int*>();
    // std::fill(bo0_map, bo0_map + DATA_SIZE, 0);
    // std::fill(bo1_map, bo1_map + DATA_SIZE, 0);
    // std::fill(bo2_map, bo2_map + DATA_SIZE, 0);

    // // Create the test data
    // int bufReference[DATA_SIZE];
    // for (int i = 0; i < DATA_SIZE; ++i) {
    //     bo0_map[i] = i;
    //     bo1_map[i] = i;
    //     bufReference[i] = bo0_map[i] + bo1_map[i]; //Generate check data for validation
    // }


    // // Synchronize buffer content with device side
    // std::cout << "synchronize input buffer data to device global memory\n";
    // boIn1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    // boIn2.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    // std::cout << "Execution of the kernel\n";
    // auto run = krnl(boIn1, boIn2, boOut);
    // run.wait();

    // // Get the output;
    // std::cout << "Get the output data from the device" << std::endl;
    // boOut.sync(XCL_BO_SYNC_BO_FROM_DEVICE);




    double tolerance = 1e-3;
    bool passed = true;
    for (int j =0; j<DATA_SIZE; j++){

        if (!(abs(bufReference[j] - output_buff[j]) < tolerance)){
            cout << "error at " << j << "expected " << bufReference[j] << "but got " << output_buff[j] << endl;
            passed = false;
        }
    }

    if (passed){
        std::cout << "TEST PASSED\n";
    } else {
        std::cout << "TEST FAILED\n";
    }
    
    return 0;
}
