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
zynqmp-common-20232:~$ ./fft_host_xrt binary_container_1.xclbin

INFO:    DATA size in words  = 4096
INFO:    DATA size in bytes  = 16384
Generating Testing Data...
Testing Data Generated
PASSED:  auto my_device = xrt::device(0)
PASSED:  auto xclbin_uuid = my_device.load_xclbin(binary_container_1.xclbin)
PASSED:  auto krnl = xrt::kernel(my_device, xclbin_uuid, "fft_top:{fft_top_1}")
INFO:    Allocate Buffers in Global Memory
PASSED:  auto bo_input = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1)) (=4))
PASSED:  auto bo_output = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2)) (=4))
PASSED:  auto bo_input_map  = bo_input.map<cmpx_data_t*>()
PASSED:  auto bo_output_map = bo_output.map<cmpx_data_t*>()
Filling Buffers with input data
INFO:    synchronize input buffer data to device global memory
PASSED:  bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE)
INFO:    Execution of the kernel

INFO:    Waiting for kernels to end...

PASSED:  run.wait()
PASSED:  bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE)


////////******* TIMING RESULTS *******////////
Timing 1D FFT Kernel with input data array of size = 2048
Time Total = 0.00025224
Load Input Time  = 1.448e-05
Excution Time = 0.0002022
Load Output Time = 4.17e-06

///////////////////////////////////////////////
Validating Output results

Error at index 1748: Expected (0.279671, -399.067), Actual (0.276428, -399.06)
                     with relative real error is 0.011594
                     with relative imag error is 1.73592e-05

Tests failed at 1 indices out of 2048

Average Relative Error in the Real Componenets Across All Indicies = 0.00688571
Average Relative Error in the Imaginary Componenets Across All Indicies = 0

```
