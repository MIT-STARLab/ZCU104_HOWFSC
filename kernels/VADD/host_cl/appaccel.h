#pragma once

#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1


//#include <CL/cl_ext_xilinx.h>
#include <vector>
#include <CL/opencl.hpp>



/**
 * Programs the device
 */
bool programXCL(const char* binaryFile, cl::Program & program,cl::Context & context,cl::CommandQueue & queue);

// buffers for exchange need to be aligned to page boundaries; otherwise CL will copy to an intermediate
inline void* alignedAllocate(const size_t num)
{
	void* ptr = NULL;	
    posix_memalign(&ptr, 4096, num );
    return ptr;
}


inline void alignedFree(void*p) 
{
    free(p);
}    


#define CHECK(error, call)  call; if (error != CL_SUCCESS) { printf("%s:%d Error %d in %s\n", __FILE__, __LINE__, error, #call); exit(EXIT_FAILURE);  }


bool listPlatforms(std::vector<cl::Platform> platforms);
bool getDevices(std::vector<cl::Device> &devices, cl_int &err, const char* vendor_name = "Xilinx");
bool is_emulation();
bool is_hw_emulation();
bool is_xpr_device(const char* device_name);