/*
 * MIT STAR Lab
 * M.Subhi Abo Rdan (msubhi_a@mit.edu)
 * Last modified on April 16, 2025
 * API for FPGA Setup/Acceleration Kernels
 */

 #include <iostream>
 #include <sys/syscall.h>
 #include <time.h>
 #include <complex>
 
 #include <xrt/xrt_device.h>
 #include <xrt/xrt_kernel.h>
 #include <xrt/xrt_bo.h>
 
 #include "fpga_kernels_api.h"
 
 
 /////// Profiling Helpers ///////
 #define TIMER(label)  timespec label; syscall(SYS_clock_gettime, CLOCK_MONOTONIC, &label)
 #define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0
 
 static const char*    STR_PASSED   = "PASSED:   ";
 static const char*    STR_FAILED   = "FAILED:   ";
 static const char*    STR_INFO     = "INFO:     ";
 static const char*    STR_USAGE    = "USAGE:    ";
 static const char*    STR_RESULTS  = "RESULTS:  ";
 
 
 /////// Acceleration Contexts ///////
 struct ZCU104_Context {
     int         dev_index;
     char*       xclbinFilename;
     xrt::device device;
     xrt::uuid   xclbin_uuid;
 
     // Kernels
     Accelerated_AngularSpectrumMethod_Context* asm_context;
     // Accelerated_2DFFT_Context* 2dfft_context;
 };
 
 struct Accelerated_AngularSpectrumMethod_Context {
     xrt::kernel      kernel;
     xrt::bo          bo_input;
     xrt::bo          bo_output;
     int              n;
 };
 
 
 /////// Methods ///////
 ZCU104_Context* create_ZCU104_context(int dev_index, char* xclbinFilename)
 {
     std::cout << STR_INFO << "FPGA SETUP CODE" << std::endl;
 
     ZCU104_Context* context = new ZCU104_Context;
     context->dev_index = dev_index;
     context->xclbinFilename = xclbinFilename;
     context->device = xrt::device(dev_index);
     context->xclbin_uuid = context->device.load_xclbin(xclbinFilename);
 
     std::cout << STR_PASSED << "FPGA SETUP CODE" << std::endl;
     return context;
 }
 
 
 
 void create_angular_spectrum_context(ZCU104_Context* ZCU104_Context, int n)
 {
     if (n!=2048) {
         std::cerr << STR_FAILED << "create_angular_spectrum_context: Matrix size not supported" << std::endl;
         return;
     }
 
     std::cout << STR_INFO << "Creating Accelerated_AngularSpectrumMethod_Context" << std::endl;
     Accelerated_AngularSpectrumMethod_Context* asm_context = new Accelerated_AngularSpectrumMethod_Context;
 
     // HLS interface
     // void angular_spectrum( bool direction, data_t distance, data_t k_2, data_t delkx, cmpx_data_t *input_mat, cmpx_data_t *output_mat);
 
     // Sizes
     asm_context->n = n;
     size_t MAT_SIZE = 2 * n * n;
     std::cout << STR_INFO << "Matrix size in words  = " << MAT_SIZE        << std::endl;
     size_t size_in_bytes = MAT_SIZE * sizeof(float); 
     std::cout << STR_INFO << "Matrix size in bytes  = " << size_in_bytes   << std::endl;
 
     // Kernel
     xrt::kernel krnl = xrt::kernel(ZCU104_Context->device, ZCU104_Context->xclbin_uuid, "angular_spectrum:{angular_spectrum_1}");
     asm_context->kernel = krnl;
     std::cout << STR_PASSED << "context->kernel    = xrt::kernel(my_device, xclbin_uuid, \"angular_spectrum:{angular_spectrum_1}\")" << std::endl;
 
     // Buffers
     asm_context->bo_input  = xrt::bo(ZCU104_Context->device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(4));            // kernel argument 4: input  matrix array
     asm_context->bo_output = xrt::bo(ZCU104_Context->device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(5));            // kernel argument 5: output matrix array
     std::cout << STR_PASSED << "context->bo_input  = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(4)) (=" << krnl.group_id(4) << "))" << std::endl;
     std::cout << STR_PASSED << "context->bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(5)) (=" << krnl.group_id(5) << "))" << std::endl;
 
     ZCU104_Context->asm_context = asm_context;
     std::cout << STR_PASSED << "Creating Accelerated_AngularSpectrumMethod_Context" << std::endl;
 }
 
 
 
 void execute_angular_spectrum_kernel(ZCU104_Context* ZCU104_Context, int n, double** wavefront, double wavelength, double pixelscale, double distance)
 {
     std::cout << STR_INFO << "Executing Angular Spectrum Kernel " << std::endl;
     if (n != ZCU104_Context->asm_context->n) { std::cerr << STR_FAILED << "execute_angular_spectrum_kernel: Matrix sizes not matching" << std::endl; return;}
 
     data_t delkx = (data_t) (2.0 * M_PI / (pixelscale * n));
     data_t k = (data_t) (2.0 * M_PI / wavelength);
     data_t k_2 = (data_t) k*k;
 
     //Map the contents of the buffer object into host memory
     auto bo_input_map   = ZCU104_Context->asm_context->bo_input.map<std::complex<data_t>*>();
     auto bo_output_map  = ZCU104_Context->asm_context->bo_output.map<std::complex<data_t>*>();
     std::cout << STR_PASSED << "auto bo_input_map  = bo_input.map<cmpx_data_t*>()"  << std::endl;
     std::cout << STR_PASSED << "auto bo_output_map = bo_output.map<cmpx_data_t*>()" << std::endl;
 
     // Fill the input buffer with data
     std::cout << STR_INFO << "synchronize input buffer data to device global memory" << std::endl;
     for (int i = 0; i < n * n; i++) {
         double real = (*wavefront)[2 * i];
         double imag = (*wavefront)[2 * i + 1];
         bo_input_map[i] = std::complex<float>(static_cast<data_t>(real), static_cast<data_t>(imag));
     }
     ZCU104_Context->asm_context->bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE);
     std::cout << STR_PASSED << "bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE)" << std::endl;
 
     // Kernel Execution
     std::cout << STR_INFO <<  "Execution of the kernel" << std::endl;
     TIMER(runTimeStart);
     auto run = ZCU104_Context->asm_context->kernel(
                 false, // is inverse_fft
                 distance,
                 k_2,
                 delkx,
                 ZCU104_Context->asm_context->bo_input,
                 ZCU104_Context->asm_context->bo_output
                 );
     run.wait();
     TIMER(runTimeEnd);
     std::cout << STR_PASSED << "run.wait()" << std::endl;
     std::cout << STR_INFO << "Kernel execution time: " << ELAPSED(runTimeEnd, runTimeStart) << " seconds" << std::endl;
 
     // Retrieving the results
     std::cout << STR_INFO << "synchronize output buffer data to host memory" << std::endl;
     ZCU104_Context->asm_context->bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
 
     for (int i = 0; i < n * n; i++) {
         (*wavefront)[2 * i]     = static_cast<double>(bo_output_map[i].real());
         (*wavefront)[2 * i + 1] = static_cast<double>(bo_output_map[i].imag());
     }
 
     std::cout << STR_PASSED << "bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE)" << std::endl;
 
     std::cout << STR_PASSED << "Executing Angular Spectrum Kernel " << std::endl;
     return;
 }
 