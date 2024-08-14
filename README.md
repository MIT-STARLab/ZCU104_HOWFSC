# High-Order WaveFront Sensing and Control (HOWFSC) on Xilinx Zynq UltraScale+ MPSoC ZCU104

Welcome to the High-Order WaveFront Sensing and Control (HOWFSC) project repository! This project demonstrates the implementation of advanced wavefront sensing and control algorithms on the Xilinx Zynq UltraScale+ MPSoC ZCU104 board. 

## Features

- **1D FFT**: Implementation and performance analysis of 1D FFT algorithms.
- **2D FFT**: Implementation and performance analysis of 2D FFT algorithms.
- **Vector Addition (VADD)**: A simple reference vector addition kernel with both HLS and software implementations.
- **Matrix-Vector Multiplication (MXV)**: Matrix-vector inner product kernel implementations.
- **MFT**: Matrix Fourier Transform (in progress)
- **Optical Modeling**: System Integration (in progress)

## Directory Structure

### [FFT SW](https://github.com/MIT-STARLab/ZCU104_HOWFSC/tree/main/FFT%20SW)
- Software implementations of 1D and 2D FFT using the Cooley–Tukey algorithm.
- Includes array and vector interfaces.
- Functions for generating random testing data.
- Scripts for runtime measurement.

### [VADD](https://github.com/MIT-STARLab/ZCU104_HOWFSC/tree/main/VADD)
- Reference implementation of a vector addition kernel.
- High-Level Synthesis (HLS) implementation.
- Software host application in C++ using Xilinx Runtime Library (XRT) with both native C++ API and OpenCL.
- Synthesis, co-simulation reports, and timing analysis.

### [1D FFT](https://github.com/MIT-STARLab/ZCU104_HOWFSC/tree/main/1D%20FFT)
- FFT IP core configurations.
- High-Level Synthesis implementation of a 1D FFT kernel.
- Software host application in C++ using Xilinx Runtime Library (XRT).
- Performance measurements and timing comparisons for FPGA vs. CPU.

### [2D FFT](https://github.com/MIT-STARLab/ZCU104_HOWFSC/tree/main/2D%20FFT)
- Block diagram of the 2D FFT design.
- High-Level Synthesis implementation of a 2D FFT kernel.
- Synthesis report and latency analysis.
- Software host application in C++ using Xilinx Runtime Library (XRT).
- Performance measurements and timing comparisons for FPGA vs. CPU.

### [MXV_INNER_PRODUCT](https://github.com/MIT-STARLab/ZCU104_HOWFSC/tree/main/MXV_INNER_PRODUCT)
- HLS and OpenCL code for matrix-vector inner product kernel.

## Resulting Performance Plots

### 1D FFT
![Timing FFT](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/1D%20FFT/timing/fft_timing_plot_fpga_mac_cpu.png)

### 2D FFT
![Timing 2D FFT](https://github.com/MIT-STARLab/ZCU104_HOWFSC/blob/main/2D%20FFT/timing/fft2d_timing_plot_fpga_mac_cpu.png)



