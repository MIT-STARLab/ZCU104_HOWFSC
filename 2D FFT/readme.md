# 2D FFT

## Block Diagram

<img width="1131" alt="Screenshot 2024-07-30 at 3 22 46 PM" src="https://github.com/user-attachments/assets/2ec11359-5173-4925-85fe-4c820f033a45">

## Hardware Deployment Results

```sh
zynqmp-common-20232:~$
zynqmp-common-20232:~$ ./fft_2d_host binary_container_1.xclbin

INFO:    Matrix size in words  = 2097152
INFO:    Matrix size in bytes  = 8388608
INFO:    Generating Testing Data...
PASSED:  Testing Data Generated Successfully
PASSED:  auto my_device = xrt::device(0)
PASSED:  auto xclbin_uuid = my_device.load_xclbin(binary_container_1.xclbin)
PASSED:  auto krnl = xrt::kernel(my_device, xclbin_uuid, "fft_2d:{fft_2d_1}")
INFO:    Allocate Buffers in Global Memory
PASSED:  auto bo_input  = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1)) (=4))
PASSED:  auto bo_temp   = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2)) (=4))
PASSED:  auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(3)) (=4))
PASSED:  auto bo_input_map  = bo_input.map<cmpx_data_t*>()
PASSED:  auto bo_output_map = bo_output.map<cmpx_data_t*>()
INFO:    Filling Buffers with input data
INFO:    synchronize input buffer data to device global memory
PASSED:  bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE)
INFO:    Execution of the kernel

INFO:    Waiting for kernels to end...

PASSED:  run.wait()
PASSED:  bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE)

************************ TIMING SUMMARY ************************
Timing 2D FFT Kernel with input matrix with size (1024, 1024)
Time Total = 0.038932
Load Input Time  = 2.548e-05
Excution Time = 0.0388528
Load Output Time = 4.62e-06
****************************************************************

INFO:   Validating Output results
Error at index 2175: Expected (-64.7979, 68219.4), Actual (-65.52, 68219.4)
relative real error is 0.0111449
relative imag error is 3.43561e-07
Error at index 6056: Expected (296957, -41.4688), Actual (296957, -42.0093)
relative real error is 1.26281e-06
relative imag error is 0.0130346
    .
    .
    .
Error at index 1047021: Expected (-91.1016, 173071), Actual (-89.2002, 173071)
relative real error is 0.0208709
relative imag error is 3.70152e-06
Error at index 1048484: Expected (-116533, 149.719), Actual (-116533, 147.833)
relative real error is 4.35767e-06
relative imag error is 0.0125952

INFO:   Tests failed at 356 entries out of 1048576

INFO:   Average Relative Error in the Real Componenet Across All entries = 4.84353e-05
INFO:   Average Relative Error in the Imaginary Componenet Across All entries = 0

```
