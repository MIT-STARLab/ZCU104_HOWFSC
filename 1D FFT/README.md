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

Timing 1D FFT Kernel with input data array of size = 2048
Time Total = 0.00024893
Load Input Time  = 1.379e-05
Excution Time = 0.00019904
Load Output Time = 4.16e-06

Validating Output results
Error at index 491: Expected (-0.567155, -0.165379), Actual (-0.564816, -0.163173)
relative real error is 0.00412473
relative imag error is 0.0133424
Error at index 492: Expected (-0.56045, -0.160922), Actual (-0.562942, -0.163371)
relative real error is 0.0044457
relative imag error is 0.0152243
Error at index 942: Expected (-0.336594, -0.0200806), Actual (-0.336335, -0.0185887)
relative real error is 0.000767119
relative imag error is 0.0742959
Error at index 1003: Expected (-0.336548, -0.0044632), Actual (-0.333357, -0.00440478)
relative real error is 0.00948119
relative imag error is 0.0130876
Error at index 1004: Expected (-0.329956, -0.00374222), Actual (-0.33352, -0.00379276)
relative real error is 0.0108014
relative imag error is 0.0135066
Error at index 1018: Expected (-0.33268, -0.00138253), Actual (-0.332867, -0.00140107)
relative real error is 0.000563654
relative imag error is 0.0134081
Error at index 1019: Expected (-0.332699, -0.00105834), Actual (-0.332879, -0.00107837)
relative real error is 0.000541048
relative imag error is 0.0189232
Error at index 1021: Expected (-0.332832, -0.000724971), Actual (-0.333005, -0.00074172)
relative real error is 0.000519699
relative imag error is 0.0231029
Error at index 1023: Expected (-0.332653, -9.62466e-05), Actual (-0.332826, -0.000108004)
relative real error is 0.000518546
relative imag error is 0.122155
Error at index 1025: Expected (-0.332829, 0.00011082), Actual (-0.332819, 0.000107884)
relative real error is 2.72209e-05
relative imag error is 0.0264892
Error at index 1106: Expected (-0.336365, 0.0187988), Actual (-0.336336, 0.0185905)
relative real error is 8.61204e-05
relative imag error is 0.0110846
Error at index 1260: Expected (-0.362332, 0.0560355), Actual (-0.361366, 0.0554777)
relative real error is 0.00266823
relative imag error is 0.00995352
Error at index 1515: Expected (-0.513167, 0.144654), Actual (-0.510968, 0.142274)
relative real error is 0.00428549
relative imag error is 0.0164489
Error at index 1516: Expected (-0.509496, 0.139986), Actual (-0.511881, 0.142569)
relative real error is 0.0046816
relative imag error is 0.0184497

Tests failed at 14 indicies out of2048


```
