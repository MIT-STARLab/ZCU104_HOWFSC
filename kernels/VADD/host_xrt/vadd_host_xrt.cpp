// Copyright © 2023 Advanced Micro Devices, Inc. All rights reserved.
// SPDX-License-Identifier: MIT

#include <sys/syscall.h>      /* Definition of SYS_* constants */
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_bo.h>

#include "vadd_include.h"



/////// Timing Helpers ///////
#define TIMER(label)  timespec label; syscall(SYS_clock_gettime, CLOCK_MONOTONIC, &label)
#define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0



/////// Error Profiling Helpers ///////
static const std::string error_message =
    "Error: Result mismatch:\n"
    "i = %f CPU result = %f Device result = %f\n";

static const char*    STR_PASSED  = "PASSED:  ";
static const char*    STR_INFO    = "INFO:    ";
static const char*    STR_USAGE   = "USAGE:   ";




typedef double data_t;

int main(int argc, char* argv[]) {

  // Check command line arguments
  if(argc < 2) {
    std::cout << STR_USAGE << argv[0] <<" <xclbin> " << std::endl;
    return EXIT_FAILURE;
  }


  // Compute the size of arrays in bytes
  size_t size_in_bytes = DATA_SIZE * sizeof(data_t);
  std::cout << STR_INFO << "DATA size in words  = " << DATA_SIZE     << std::endl;
  std::cout << STR_INFO << "DATA size in bytes  = " << size_in_bytes << std::endl;

  // Loading the xclbin
  char* xclbinFilename = argv[1];

  //open the device and return its handle
  unsigned int dev_index = 0;
  auto my_device = xrt::device(dev_index);
  std::cout << STR_PASSED << "auto my_device = xrt::device(" << dev_index << ")" << std::endl;


  //load the xclbin file from disk. This will return the UUID of the xclbin.
  auto xclbin_uuid = my_device.load_xclbin(xclbinFilename);
  std::cout << STR_PASSED << "auto xclbin_uuid = my_device.load_xclbin(" << xclbinFilename << ")" << std::endl;


  // Starting the (HLS) PL kernel
  auto krnl = xrt::kernel(my_device, xclbin_uuid, "krnl_vadd:{krnl_vadd_1}");
  std::cout << STR_PASSED << "auto krnl = xrt::kernel(my_device, xclbin_uuid, \"krnl_vadd:{krnl_vadd_1}\")" << std::endl;


  std::cout << STR_INFO << "Allocate Buffers in Global Memory" << std::endl;
  //The buffer type is default buffer object (bo) with host buffer and device buffer. The host buffer is allocated and managed by XRT.
  auto bo0 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(0));  //input 0
  auto bo1 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1));  //input 1
  auto bo2 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2));  //output
  std::cout << STR_PASSED << "auto bo0 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(0) (=" << krnl.group_id(0) << "))" << std::endl;
  std::cout << STR_PASSED << "auto bo1 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1) (=" << krnl.group_id(1) << "))" << std::endl;
  std::cout << STR_PASSED << "auto bo2 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2) (=" << krnl.group_id(2) << "))" << std::endl;

  //Map the contents of the buffer object into host memory
  auto bo0_map = bo0.map<data_t*>();
  auto bo1_map = bo1.map<data_t*>();
  auto bo2_map = bo2.map<data_t*>();
  std::fill(bo0_map, bo0_map + DATA_SIZE, 0);
  std::fill(bo1_map, bo1_map + DATA_SIZE, 0);
  std::fill(bo2_map, bo2_map + DATA_SIZE, 0);
  std::cout << STR_PASSED << "auto bo0_map = bo0.map<data_t*>()" << std::endl;
  std::cout << STR_PASSED << "auto bo1_map = bo1.map<data_t*>()" << std::endl;
  std::cout << STR_PASSED << "auto bo2_map = bo2_map<data_t*>()" << std::endl;

  //setting input data
  std::cout << STR_INFO << "Creating random input data" << std::endl;
  FILE *fin1=fopen("vector_inputs.txt", "w");
  srand(time(NULL));   // Initialization, should only be called once.
  for(int i = 0 ; i< DATA_SIZE; i++) {
      data_t r1 = (data_t) rand(); // Returns a pseudo-random integer between 0 and RAND_MAX.
      data_t r2 = (data_t) rand(); // Returns a pseudo-random integer between 0 and RAND_MAX.
      bo0_map[i] = r1;
      bo1_map[i] = r2;
      fprintf(fin1, "%10f\t%10f\n", r1, r2);
    }
  fclose(fin1);



  // Synchronize buffer content with device side
  std::cout << STR_INFO << "synchronize input buffer data to device global memory" << std::endl;

  TIMER(startTime);

  bo0.sync(XCL_BO_SYNC_BO_TO_DEVICE);
  std::cout << STR_PASSED << "bo0.sync(XCL_BO_SYNC_BO_TO_DEVICE)" << std::endl;
  bo1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
  std::cout << STR_PASSED << "bo1.sync(XCL_BO_SYNC_BO_TO_DEVICE)" << std::endl;

  TIMER(loadTime);

  std::cout << STR_INFO <<  "Execution of the kernel" << std::endl;
  auto run = krnl(bo0, bo1, bo2, DATA_SIZE);
  // Waiting for all kernels to end
  std::cout << std::endl << STR_INFO << "Waiting for kernels to end..." << std::endl << std::endl;
  run.wait();
  std::cout << STR_PASSED << "run.wait()" << std::endl;

  TIMER(runTime);


  // Retrieving the results
  bo2.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
  std::cout << STR_PASSED << "bo2.sync(XCL_BO_SYNC_BO_FROM_DEVICE)" << std::endl;

  TIMER(doneTime);

  printf("Time total %f load %f run %f readout %f\n", ELAPSED(doneTime,startTime),ELAPSED(loadTime,startTime),ELAPSED(runTime,loadTime),ELAPSED(doneTime,runTime));


  // Verify the result
  std::cout << STR_INFO <<  "Checking the results" << std::endl;
  FILE *fin2=fopen("vector_out.txt", "w");
  int match = 0;
  for (int i = 0; i < DATA_SIZE; i++)
  {
      data_t host_result = bo0_map[i] + bo1_map[i];
      fprintf(fin2, "%10f\t%10f\n", host_result, bo2_map[i]);
      if (bo2_map[i] != host_result)
	    {
    	  printf(error_message.c_str(), i, host_result, bo2_map[i]);
    	  match = 1;
    	  break;
	    }
  }
  fclose(fin2);

  std::cout << "TEST " << (match ? "FAILED" : "PASSED") << std::endl;
  return (match ? EXIT_FAILURE : EXIT_SUCCESS);
}
