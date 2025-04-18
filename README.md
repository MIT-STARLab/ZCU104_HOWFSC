# High-Order WaveFront Sensing and Control (HOWFSC) on Xilinx ZCU104

Welcome to the High-Order WaveFront Sensing and Control (HOWFSC) project repository!  
This project demonstrates the implementation and performance evaluation of advanced wavefront sensing and control algorithms on the Xilinx Zynq UltraScale+ MPSoC ZCU104 board using High-Level Synthesis (HLS).


## 🛰️ Overview

This project focuses on accelerating computational bottlenecks in Electric Field Conjugation (EFC) and Jacobian Estimation algorithms, which are critical in high-order wavefront sensing and control systems. By leveraging HLS, we implement and optimize key computational kernels to run efficiently on the ZCU104 platform.



## 📁 Repository Structure

```text
ZCU104_HOWFSC/
│
├── kernels/                     # All HLS computational kernels including HLS code, host XRT code, test benches, sw implementations, design block diagrams, performance analysis
│   ├── 1D FFT/                          # 1D FFT kernel + reports
│   ├── 2D FFT/                          # 2D FFT kernel + reports
│   ├── angular_spectrum_propagation/    # Angular Spectrum Method
│   ├── MFT/                             # Matrix Fourier Transform  (in progress)
│   ├── MXV/                             # Matrix-Vector Inner Product
│   ├── QR_Givens/                       # Given's QR decomposition
│   └── VADD/                            # Simple vector add baseline
│
├── accelerated_howfsc           # Full Embedded Implementation of EFC and Optical Modeling with Host Server (Kernels Integration in Progress)
└── README.md                            # This file
```


## 👏 Acknowledgments
This project is part of the NASA APRA-funded HOWFSC effort at the MIT STAR Lab
