/*
 * MIT STAR Lab
 * M.Subhi Abo Rdan (msubhi_a@mit.edu)
 * Last modified on April 16, 2025
 * API for FPGA Setup/Acceleration Kernels
 */

#ifndef FPGA_KERNELS_API_H
#define FPGA_KERNELS_API_H

#ifdef __cplusplus
extern "C" { // Extern as a C API
#endif


typedef float data_t;                                               // data type for the kernel: float only
typedef struct ZCU104_Context ZCU104_Context;                       // defined in fpga_kernels_api.cpp
typedef struct Accelerated_AngularSpectrumMethod_Context Accelerated_AngularSpectrumMethod_Context;     // defined in fpga_kernels_api.cpp


/**
 * @brief Constructor:   ZCU104 context for FPGA acceleration
 * @param dev_index      The index of the FPGA device to use (default is 0)
 * @param xclbinFilename The path to the xclbin file (default is "howfsc_kernels.xclbin")
 */
ZCU104_Context* create_ZCU104_context(int dev_index , char* xclbinFilename);


/***************       Methods     ******************/

/**
 * @brief Create a context for the Angular Spectrum Method kernel
 * @param ZCU104_Context The ZCU104 context created by create_ZCU104_context
 * @param n              The size of the matrix
 */
void create_angular_spectrum_context(ZCU104_Context* ZCU104_Context, int n);


/**
 * @brief Execute the Angular Spectrum Method kernel
 * @param context        The context created by create_angular_spectrum_context
 * @param n              The size of the matrix (must match context->n)
 * @param wavefront      Pointer to Pointer to the input wavefront 2D array of complex numbers
 * @param wavelength     The wavelength of the light
 * @param pixelscale     The pixel scale of the wavefront
 * @param distance       The distance to propagate the wavefront
 */
void execute_angular_spectrum_kernel(ZCU104_Context* ZCU104_Context, int n, double** wavefront, double wavelength, double pixelscale, double distance);


#ifdef __cplusplus
}
#endif

#endif //FPGA_KERNELS_API_H