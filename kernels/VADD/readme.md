# VADD HLS

## Optimizations

### `#pragma HLS INTERFACE m_axi`
- The `m_axi` interface in HLS provides independent read and write channels, allowing for simultaneous read and write operations.
- This setup enables reading input data for the next cycle (i+1) while writing output data of the current cycle (i), enhancing pipelining efficiency.
- **gmem0** is used for both **in2** and **out**.
- Placing both **in1** and **in2** on the same interface increases latency due to the single read operation limitation per interface.

By parallelizing the reading of inputs, we achieve improved performance at the cost of additional resources.

An alternative optimization could involve using the AXI stream (`axis`) interface, which is beneficial for sequential data reading. However, this requires the kernel to interact with another streaming component. In our scenario, the kernel interacts with the host application, making the use of the stream interface unsuitable. Therefore, we stick with the memory-mapped interface.

### `#pragma HLS DATAFLOW`
The `dataflow` pragma instructs the compiler to run the following three APIs in parallel, allowing `in1` and `in2` to load simultaneously since they are independent.

### `#pragma HLS PIPELINE`
The `pipeline` pragma allows starting the next iteration in the loop before finishing the current iteration.

## Timing Analysis

![Timing Analysis](https://github.com/user-attachments/assets/7e031e93-0985-465b-b9b7-8a281df56e82)

For all the loops, we can observe that the total latency for the loop is approximately `DATA_SIZE + iteration latency`.

For example:
- The load input loop has an iteration latency of 3 cycles and an interval of 1 cycle. This means that when we start iteration `i`, we request the loading element at index `i` from memory and then start the next iteration in the next cycle.
- Each iteration takes 3 cycles, resulting in a total latency for the loop of approximately `DATA_SIZE + latency = 4172`.

Notice how the total latency of the kernel is almost the same as the latency of each component, reflecting how the kernels are working in parallel using the `#pragma DATAFLOW`. If we removed this pragma, we would notice the latency add up.

![Latency Comparison](https://github.com/user-attachments/assets/7dacc473-2702-47e5-8f62-01859ad999db)

## Synthesis Report Warnings

- **WARNING: [HLS 200-805] An internal stream 'out_stream' with default size can result in deadlock. Please consider resizing the stream using the directive `set_directive_stream` or the `HLS stream` pragma.** 
  - While our kernel doesn't cause deadlocks, it may reduce throughput slightly by stalling `input_load` because the size of the FIFO in the streams is only 2 and the `load_input` loop latency is 3, but `compute` is 6. The compute is slower. We can resize the FIFO with resources tradeoff.

- **WARNING: [RTGEN 206-101] Design contains AXI ports. Reset is fixed to synchronous and active low.** 
  - This is acceptable as we are using AXI.

- **WARNING: [HLS 200-656] Deadlocks can occur since process `load_input` is instantiated in a dataflow region with `ap_ctrl_none` or without start propagation and contains an auto-rewind pipeline.**
  - This warning disappears when we stop the pipeline inside the loops. However, it is not causing any problems.

Other warnings do not seem important.


## Hardware Deployment Results
```sh
zynqmp-common-20232:~$ ls
BOOT.BIN  Image  binary_container_1.xclbin  boot.scr  vadd_host_xrt

zynqmp-common-20232:~$ ./vadd_host_xrt binary_container_1.xclbin
INFO:    DATA size in words  = 4096
INFO:    DATA size in bytes  = 32768
PASSED:  auto my_device = xrt::device(0)
PASSED:  auto xclbin_uuid = my_device.load_xclbin(binary_container_1.xclbin)
PASSED:  auto krnl = xrt::kernel(my_device, xclbin_uuid, "krnl_vadd:{krnl_vadd_1}")
INFO:    Allocate Buffers in Global Memory
PASSED:  auto bo0 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(0) (=4))
PASSED:  auto bo1 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(1) (=4))
PASSED:  auto bo2 = xrt::bo(my_device, size_in_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(2) (=4))
PASSED:  auto bo0_map = bo0.map<data_t*>()
PASSED:  auto bo1_map = bo1.map<data_t*>()
PASSED:  auto bo2_map = bo2_map<data_t*>()
INFO:    Creating random input data
INFO:    synchronize input buffer data to device global memory
PASSED:  bo0.sync(XCL_BO_SYNC_BO_TO_DEVICE)
PASSED:  bo1.sync(XCL_BO_SYNC_BO_TO_DEVICE)
INFO:    Execution of the kernel

INFO:    Waiting for kernels to end...

PASSED:  run.wait()
PASSED:  bo2.sync(XCL_BO_SYNC_BO_FROM_DEVICE)
Time total 0.000301 load 0.000039 run 0.000248 readout 0.000014
INFO:    Checking the results
TEST PASSED

zynqmp-common-20232:~$
```

