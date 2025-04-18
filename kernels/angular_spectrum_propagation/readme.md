# Angular Spectrum Propagation

## Design Block Diagram
![Block Diagram](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/kernels/angular_spectrum_propagation/angular_spectrum_kernel_block_diagram.png)


## Timing
![Timing](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/kernels/angular_spectrum_propagation/ang_spec_timing_plot_fpga_mac_cpu.png)
~ 465 times faster than CPU




## Hardware Deployment Results
![Block Diagram](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/kernels/angular_spectrum_propagation/validation/sig1.png)



```sh
zynqmp-common-20232:~$ ./angular_spectrum_host binary_container_1.xclbin 5 0                                                                                                                             
PASSED:   auto my_device = xrt::device(0)
PASSED:   auto xclbin_uuid = my_device.load_xclbin(binary_container_1.xclbin)
PASSED:   auto krnl = xrt::kernel(my_device, xclbin_uuid, "angular_spectrum:{angular_spectrum_1}")
INFO:     Allocate Buffers in Global Memory
INFO:     Matrix size in words  = 2097152
INFO:     Matrix size in bytes  = 8388608
PASSED:   auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(4)) (=4))
PASSED:   auto bo_ouput = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(5)) (=4))
PASSED:   auto bo_input_map  = bo_input.map<cmpx_data_t*>()
PASSED:   auto bo_output_map = bo_output.map<cmpx_data_t*>()
INFO:     Generating Testing Data...
PASSED:   Testing Data Generated
INFO:     Filling Argument Buffers with input data
INFO:     Running on CPU First
PASSED:   CPU Done
INFO:     synchronize input buffer data to device global memory
PASSED:   bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE)
INFO:     Execution of the kernel
PASSED:   run.wait()
PASSED:   bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE)
Round 0 Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.277521
Total FPGA Time       =   0.024250
->Input Load Time     =   0.000036
->Excution Time       =   0.024164
->Output Load Time    =   0.000004

INFO:     RUN KERNEL NUMBER OF TIMES =  5
INFO:     STARTING ROUND : 0
Round 1 Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.277707
Total FPGA Time       =   0.000005
->Input Load Time     =   0.000021
->Excution Time       =   0.023965
->Output Load Time    =   0.023992

INFO:     STARTING ROUND : 1
Round 2 Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.278197
Total FPGA Time       =   0.000005
->Input Load Time     =   0.000020
->Excution Time       =   0.023950
->Output Load Time    =   0.023975

INFO:     STARTING ROUND : 2
Round 3 Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.277986
Total FPGA Time       =   0.000005
->Input Load Time     =   0.000021
->Excution Time       =   0.023945
->Output Load Time    =   0.023971

INFO:     STARTING ROUND : 3
Round 4 Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.280425
Total FPGA Time       =   0.000005
->Input Load Time     =   0.000021
->Excution Time       =   0.023981
->Output Load Time    =   0.024006

INFO:     STARTING ROUND : 4
Round 5 Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   11.277357
Total FPGA Time       =   0.000006
->Input Load Time     =   0.000021
->Excution Time       =   0.023955
->Output Load Time    =   0.023981


**************************      Timing Summary   ****************************
Kernel name                               = Angular Spectrum Propagation
Input Data Size                           = (1024, 1024)
Number of Run Rounds                      = 5
FPGA Geometric Mean Run Time (With Sync)  = 0.023985
FPGA Geometric Mean Run Time (No   Sync)  = 0.023959
CPU  Geometric Mean Run Time              = 11.278334
FPGA Arithmetic Mean Run Time (With Sync)  = 0.023985
FPGA Arithmetic Mean Run Time (No   Sync)  = 0.023959
CPU  Arithmetic Mean Run Time              = 11.278335
*****************************************************************************
```
