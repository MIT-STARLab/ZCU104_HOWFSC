# Angular Spectrum Propagation

## Design Block Diagram
![Block Diagram](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/angular_spectrum_propagation/angular_spectrum_kernel_block_diagram.png)


## Hardware Deployment Results
![Block Diagram](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/angular_spectrum_propagation/deployment_results.png)

~ 465 times faster than CPU

```sh
zynqmp-common-20232:~$ ./angular_spectrum_host binary_container_1.xclbin 5 0                                                                                                                             
PASSED:   auto my_device = xrt::device(0)
PASSED:   auto xclbin_uuid = my_device.load_xclbin(binary_container_1.xclbin)
PASSED:   auto krnl = xrt::kernel(my_device, xclbin_uuid, "angular_spectrum:{angular_spectrum_1}")
INFO:     Allocate Buffers in Global Memory
INFO:     Matrix size in words  = 2097152
INFO:     Matrix size in bytes  = 8388608
PASSED:   auto bo_kxy   = xrt::bo(my_device, MAT_ROWS * sizeof(data_t), XCL_BO_FLAGS_NONE, krnl.group_id(4)) (=65535))
PASSED:   auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(5)) (=4))
PASSED:   auto bo_ouput = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(6)) (=4))
PASSED:   auto bo_kxy_map    = bo_kxy.map<data_t*>()
PASSED:   auto bo_input_map  = bo_input.map<cmpx_data_t*>()
PASSED:   auto bo_output_map = bo_output.map<cmpx_data_t*>()
INFO:     Generating Testing Data...
PASSED:   Testing Data Generated
INFO:     Filling Argument Buffers with input data
INFO:     Running on CPU First
PASSED:   CPU Done
INFO:     synchronize input buffer data to device global memory
PASSED:   bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE)
PASSED:   bo_kxy.sync(XCL_BO_SYNC_BO_TO_DEVICE)
INFO:     Execution of the kernel
PASSED:   run.wait()
PASSED:   bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE)

Round Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.277367
Total FPGA Time       =   0.024287
->Input Load Time     =   0.000041
->Excution Time       =   0.024186
->Output Load Time    =   0.000004

INFO:     RUN KERNEL NUMBER OF TIMES =  5
INFO:     STARTING ROUND : 0
Round Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.302866
Total FPGA Time       =   0.000006
->Input Load Time     =   0.000024
->Excution Time       =   0.023953
->Output Load Time    =   0.023983
INFO:     Start Validation

INFO:     STARTING ROUND : 1
Round Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.302166
Total FPGA Time       =   0.000006
->Input Load Time     =   0.000025
->Excution Time       =   0.023959
->Output Load Time    =   0.023990
INFO:     Start Validation

INFO:     STARTING ROUND : 2
Round Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.302623
Total FPGA Time       =   0.000006
->Input Load Time     =   0.000025
->Excution Time       =   0.023950
->Output Load Time    =   0.023981
INFO:     Start Validation

INFO:     STARTING ROUND : 3
Round Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.302531
Total FPGA Time       =   0.000005
->Input Load Time     =   0.000026
->Excution Time       =   0.023964
->Output Load Time    =   0.023995
INFO:     Start Validation

INFO:     STARTING ROUND : 4
Round Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.302573
Total FPGA Time       =   0.000005
->Input Load Time     =   0.000024
->Excution Time       =   0.023940
->Output Load Time    =   0.023969

**************************      Timing Summary   ****************************
Kernel name                               = Angular Spectrum Propagation
Input Data Size                           = (1024, 1024)
Number of Run Rounds                      = 5
FPGA Geometric Mean Run Time (With Sync)  = 0.023984
FPGA Geometric Mean Run Time (No   Sync)  = 0.023953
CPU  Geometric Mean Run Time              = 11.302552
FPGA Arithmetic Mean Run Time (With Sync)  = 0.023984
FPGA Arithmetic Mean Run Time (No   Sync)  = 0.023953
CPU  Arithmetic Mean Run Time              = 11.302552
*****************************************************************************
```