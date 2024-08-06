# 1D FFT Kernel

## FFT IP Core Static Configuration Parameters
This section describes the configuration parameters used for the HLS FFT implementation.

| **Parameter**                  | **Description**                                     | **Value**    |
|--------------------------------|-----------------------------------------------------|--------------|
| `DATA_SIZE`                    | Size of the data array                             | 2048         |
| `FFT_INPUT_WIDTH`              | Bit-width of the input data (e.g., size of `float`) | 32           |
| `FFT_OUTPUT_WIDTH`             | Bit-width of the output data                       | 32           |
| `FFT_NFFT_MAX`                 | Maximum FFT length (log2 of `DATA_SIZE`)           | 11           |
| `HAS_NFFT`                     | Runtime configurable length flag                   | false        |
| `FFT_PHASE_FACTOR_WIDTH`       | Precision of the internal phase factor              | 24           |
| `ORDERING_OPT`                 | Output data order (natural or bit-reversed)        | `natural_order` |
| `FFT_STAGES_BLOCK_RAM`         | Number of block RAM stages used                    | 2            |

For more notes on the general static parameters, go to the readme file in the [**/hls**](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/1D%20FFT/hls/readme.md) folder


# HLS Timing and Constraints

### Synthesis Report
![Screenshot 2024-07-29 at 3 29 08 PM](https://github.com/user-attachments/assets/683052db-1eab-4113-8d65-cb5514644c64)

### Co-Simulation Report
![Screenshot 2024-07-29 at 3 45 10 PM](https://github.com/user-attachments/assets/d1485a1d-5df8-4e16-93e7-ed3a7d1101b6)

## Co-Simulation Results Graph
The data in the test bench was saved to files to be plotted using the [python script](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/1D%20FFT/hls/plot_sim_results.py).
![Screenshot 2024-07-29 at 2 58 48 PM](https://github.com/user-attachments/assets/157e88f5-b701-4e17-94e5-edf9dc1ebf04)

## Hardware Deployment Results
```sh
zynqmp-common-20232:~$
zynqmp-common-20232:~$ ./fft_host binary_container_1.xclbin 10

PASSED:    auto my_device = xrt::device(0)
PASSED:    auto xclbin_uuid = my_device.load_xclbin(binary_container_1.xclbin)
PASSED:    auto krnl = xrt::kernel(my_device, xclbin_uuid, "fft_top:{fft_top_1}")

INFO:  Allocate Buffers in Global Memory
INFO:  DATA size in words  = 4096
INFO:  DATA size in bytes  = 16384
PASSED:    auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1)) (=4))
PASSED:    auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2)) (=4))
PASSED:    auto bo_input_map  = bo_input.map<cmpx_data_t*>()
PASSED:    auto bo_output_map = bo_output.map<cmpx_data_t*>()

INFO:  Generating Testing Data...
PASSED:    Testing Data Generated
INFO:  Filling Argument Buffers with input data

INFO:  Running on CPU First
PASSED:    CPU Done

INFO:  synchronize input buffer data to device global memory
PASSED:    bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE)
INFO:  Execution of the kernel

INFO:  Waiting for kernels to end...

PASSED:    run.wait()
PASSED:    bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE)

INFO:  Start Validation
Error   Results Mismatch:
        Index i = 1748; CPU Results (0.279671, 399.066528); Device Result (0.276428, 399.059631)
        Relative error in real component = 0.011594
        Relative error in imag component = 0.000017

RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 4.57271e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 2.69106e-05
FAILED: Validation; test failed at 1 indicies out of 2048
**************************   Timing Summary **************************
Data Size           =   2048
Total CPU Time      =   0.003783
Total FPGA Time     =   0.000267
->Input Load Time   =   0.000014
->Excution Time     =   0.000217
->Output Load Time  =   0.000004
***********************************************************************

INFO:  RUN KERNEL NUMBER OF TIMES =  10

INFO:  STATING ROUND : 0
INFO:  Start Validation
Error:  Results Mismatch:
        Index i = 147; CPU Results (17556.531250, 1.573242); Device Result (17556.517578, 1.600830)
        Relative error in real component = 0.000001
        Relative error in imag component = 0.017536
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 1.00418e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 2.23513e-05
FAILED: Validation; test failed at 1indicies out of 2048
**************************   Timing Summary **************************
Data Size           =   2048
Total CPU Time      =   0.003792
Total FPGA Time     =   0.000108
->Input Load Time   =   0.000005
->Excution Time     =   0.000098
->Output Load Time  =   0.000004
***********************************************************************

INFO:  STATING ROUND : 1
INFO:  Start Validation
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 2.08419e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 1.50779e-05
PASSED:    Validation; data at all indicies match
**************************   Timing Summary **************************
Data Size           =   2048
Total CPU Time      =   0.003772
Total FPGA Time     =   0.000109
->Input Load Time   =   0.000004
->Excution Time     =   0.000100
->Output Load Time  =   0.000004
***********************************************************************

INFO:  STATING ROUND : 2
INFO:  Start Validation
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 1.25822e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 1.21944e-05
PASSED:    Validation; data at all indicies match
**************************   Timing Summary **************************
Data Size           =   2048
Total CPU Time      =   0.003781
Total FPGA Time     =   0.000107
->Input Load Time   =   0.000004
->Excution Time     =   0.000097
->Output Load Time  =   0.000004
***********************************************************************

INFO:  STATING ROUND : 3
INFO:  Start Validation
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 1.08558e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 1.13566e-05
PASSED:    Validation; data at all indicies match
**************************   Timing Summary **************************
Data Size           =   2048
Total CPU Time      =   0.003796
Total FPGA Time     =   0.000107
->Input Load Time   =   0.000004
->Excution Time     =   0.000097
->Output Load Time  =   0.000004
***********************************************************************

INFO:  STATING ROUND : 4
INFO:  Start Validation
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 8.31326e-06
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 2.04797e-05
PASSED:    Validation; data at all indicies match
**************************   Timing Summary **************************
Data Size           =   2048
Total CPU Time      =   0.003782
Total FPGA Time     =   0.000106
->Input Load Time   =   0.000004
->Excution Time     =   0.000097
->Output Load Time  =   0.000004
***********************************************************************

INFO:  STATING ROUND : 5
INFO:  Start Validation
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 1.25544e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 1.0356e-05
PASSED:    Validation; data at all indicies match
**************************   Timing Summary **************************
Data Size           =   2048
Total CPU Time      =   0.003780
Total FPGA Time     =   0.000107
->Input Load Time   =   0.000004
->Excution Time     =   0.000097
->Output Load Time  =   0.000004
***********************************************************************

INFO:  STATING ROUND : 6
INFO:  Start Validation
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 1.57015e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 1.07352e-05
PASSED:    Validation; data at all indicies match
**************************   Timing Summary **************************
Data Size           =   2048
Total CPU Time      =   0.003783
Total FPGA Time     =   0.000107
->Input Load Time   =   0.000004
->Excution Time     =   0.000097
->Output Load Time  =   0.000004
***********************************************************************

INFO:  STATING ROUND : 7
INFO:  Start Validation
Error:  Results Mismatch:
        Index i = 891; CPU Results (-5198.112793, 5.483093); Device Result (-5198.113770, 5.416473)
        Relative error in real component = 0.000000
        Relative error in imag component = 0.012150
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 9.40051e-06
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 1.94657e-05
FAILED: Validation; test failed at 1indicies out of 2048

**************************   Timing Summary **************************
Data Size           =   2048
Total CPU Time      =   0.003779
Total FPGA Time     =   0.000106
->Input Load Time   =   0.000004
->Excution Time     =   0.000096
->Output Load Time  =   0.000004
***********************************************************************

INFO:  STATING ROUND : 8
INFO:  Start Validation
Error:  Results Mismatch:
        Index i = 1944; CPU Results (0.053711, 7781.453613); Device Result (0.087036, 7781.503418)
        Relative error in real component = 0.620455
        Relative error in imag component = 0.000006
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 0.000327228
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 1.10576e-05
FAILED: Validation; test failed at 1indicies out of 2048

**************************   Timing Summary **************************
Data Size           =   2048
Total CPU Time      =   0.003780
Total FPGA Time     =   0.000107
->Input Load Time   =   0.000004
->Excution Time     =   0.000098
->Output Load Time  =   0.000004
***********************************************************************

INFO:  STATING ROUND : 9
INFO:  Start Validation
RESULTS: Average Relative Error in the Real Componenet Across All Indicies = 1.32928e-05
RESULTS: Average Relative Error in the Imag Componenet Across All Indicies = 1.19009e-05
PASSED:    Validation; data at all indicies match

**************************   Timing Summary **************************
Data Size           =   2048
Total CPU Time      =   0.003779
Total FPGA Time     =   0.000106
->Input Load Time   =   0.000004
->Excution Time     =   0.000097
->Output Load Time  =   0.000004
***********************************************************************

PASSED:    RUN KERNEL NUMBER OF TIMES =  10


****************   Timing Summary **************************
************************************************************
RESULTS: FPGA Time average  = 0.000107018
RESULTS: CPU  Time average  = 0.00378236
************************************************************
************************************************************

```
