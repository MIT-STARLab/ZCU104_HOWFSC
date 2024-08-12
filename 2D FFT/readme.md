# 2D FFT

## Block Diagram
![Block Diagram](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/2D%20FFT/block_diagram_v1.png)

## Performance and Resources (Synthesis Report)
![Synthesis Report](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/2D%20FFT/syn_rpt.png)



## Timing Summary for 2D FFT

The following table provides a detailed summary of the timing results for the 2D FFT kernel run on FPGA and CPU. 
The measurements are taken for different input data sizes, and both with and without synchronization for FPGA.
Testing data is generated randomly ~Unif(100, 1000)

| FFT Size | Latency (cycles)            | FPGA Average Run Time (With Sync) | FPGA Average Run Time (No Sync) | CPU Average Run Time |
|----------|-----------------------------|-----------------------------------|---------------------------------|-----------------------|
| 128      | 31,119                      | 0.000311 seconds     | 0.000301 seconds                | 0.041216 seconds      |
| 256      | 112,079                     | 0.001005 seconds     | 0.186729 seconds                | 0.186729 seconds      |
| 512      | 423,503                     | 0.003727 seconds     | 1.144043 seconds                | 1.144043 seconds      |
| 1024     | 1,636,431                   | 0.015062 seconds     | 5.142044 seconds                | 5.142044 seconds      |
| 2048     | 6,429,775                   | 0.065410 seconds     | 23.538236 seconds               | 23.538236 seconds     |


![Timing](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/2D%20FFT/timing/fft2d_timing_plot_fpga_mac_cpu.png)


**Notes:**
- **Latency** is provided in cycles and is calculated from Vitis synthesis report. Convert using * (10 ns/cycle) * (10^-9 ns/s).
- **FPGA Average Run Time (With Sync)** includes synchronization overhead.
- **FPGA Average Run Time (No Sync)** excludes synchronization overhead.
- **CPU Average Run Time** is the run time on the ZCU 104 CPU for the same FFT size.
- **Mac Timing Data** are in the [FFT SW](https://github.com/MIT-STARLab/ZCU104_HOWFSC/tree/main/FFT%20SW)/timing Folder




## Hardware Deployment Results

```sh
zynqmp-common-20232:~$ 
zynqmp-common-20232:~$ ./fft_2d_host binary_container_1.xclbin  

PASSED:    auto my_device = xrt::device(0)
PASSED:    auto xclbin_uuid = my_device.load_xclbin(binary_container_1.xclbin)
PASSED:    auto krnl = xrt::kernel(my_device, xclbin_uuid, "fft_2d:{fft_2d_1}")
INFO:  Allocate Buffers in Global Memory
INFO:  Matrix size in words  = 2097152
INFO:  Matrix size in bytes  = 8388608
PASSED:    auto bo_input  = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1)) (=4))
PASSED:    auto bo_temp   = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2)) (=4))
PASSED:    auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(3)) (=4))
PASSED:    auto bo_input_map  = bo_input.map<cmpx_data_t*>()
PASSED:    auto bo_output_map = bo_output.map<cmpx_data_t*>()
INFO:  Generating Testing Data...
PASSED:    Testing Data Generated
INFO:  Filling Argument Buffers with input data
Filling Buffers with input data
INFO:  Running on CPU First
PASSED:    CPU Done
INFO:  synchronize input buffer data to device global memory
PASSED:    bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE)
INFO:  Execution of the kernel

INFO:  Waiting for kernels to end...

PASSED:    run.wait()
PASSED:    bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE)

INFO:  Start Validation

Error:  Results Mismatch:
        Index i = 14925; CPU Results (-340990.093750, 19.281250); Device Result (-340990.343750, 19.581055)
        Relative error in real component = 0.000001
        Relative error in imag component = 0.015549
...
Error:  Results Mismatch:
        Index i = 1047184; CPU Results (6.187500, 13021.396484); Device Result (4.573242, 13022.248047)
        Relative error in real component = 0.260890
        Relative error in imag component = 0.000065

RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 2.70611e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 2.98263e-05
FAILED: Validation; test failed at 213 indicies out of 1048576
2D FFT Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   5.193193
Total FPGA Time       =   0.015350
->Input Load Time     =   0.000033
->Excution Time       =   0.015262
->Output Load Time    =   0.000005



*********************************************
INFO:  RUN KERNEL NUMBER OF TIMES =  10
*********************************************

INFO:  STATING ROUND : 0
INFO:  Start Validation
Error:  Results Mismatch:
        Index i = 291; CPU Results (2.328125, 394607.500000); Device Result (1.496582, 394607.781250)
        Relative error in real component = 0.357173
        Relative error in imag component = 0.000001
.
.
.
Error:  Results Mismatch:
        Index i = 1047584; CPU Results (169052.171875, 27.906250); Device Result (169050.265625, 27.387207)
        Relative error in real component = 0.000011
        Relative error in imag component = 0.018600
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 3.2744e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 2.66638e-05
FAILED: Validation; test failed at 164 indicies out of 1048576
2D FFT Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   5.192590
Total FPGA Time       =   0.015082
->Input Load Time     =   0.000019
->Excution Time       =   0.015056
->Output Load Time    =   0.000005

INFO:  STATING ROUND : 1
INFO:  Start Validation
Error:  Results Mismatch:
        Index i = 3412; CPU Results (14.884766, 138721.921875); Device Result (15.179688, 138722.171875)
        Relative error in real component = 0.019814
        Relative error in imag component = 0.000002
.
.
.
Error:  Results Mismatch:
        Index i = 1046236; CPU Results (59.882812, -4783.828125); Device Result (57.494141, -4785.624023)
        Relative error in real component = 0.039889
        Relative error in imag component = 0.000375
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 5.85683e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 3.34888e-05
FAILED: Validation; test failed at 178indicies out of 1048576
2D FFT Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   5.192721
Total FPGA Time       =   0.015084
->Input Load Time     =   0.000019
->Excution Time       =   0.015057
->Output Load Time    =   0.000005
.
.
.
.
.
.

INFO:  STATING ROUND : 9
INFO:  Start Validation
Error:  Results Mismatch:
        Index i = 9394; CPU Results (7.539062, 158312.906250); Device Result (7.638184, 158312.593750)
        Relative error in real component = 0.013148
        Relative error in imag component = 0.000002
...
Error:  Results Mismatch:
        Index i = 1028304; CPU Results (-173891.093750, 32.265625); Device Result (-173889.843750, 34.137207)
        Relative error in real component = 0.000007
        Relative error in imag component = 0.058005
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 2.10183e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 2.47206e-05
FAILED: Validation; test failed at 165indicies out of 1048576
2D FFT Run Time:
-------------------------
Matrix Dimension Size =   (1024, 1024)
Total CPU Time        =   5.207527
Total FPGA Time       =   0.015087
->Input Load Time     =   0.000019
->Excution Time       =   0.015060
->Output Load Time    =   0.000005




**************************      Timing Summary   ****************************
Kernel name                        = 2D FFT
Input Data Size                    = (1024, 1024)
Number of Run Rounds               = 10
FPGA Average Run Time (With Sync)  = 0.015082
FPGA Average Run Time (No   Sync)  = 0.015056
CPU  Average Run Time              = 5.200041
*****************************************************************************

```
