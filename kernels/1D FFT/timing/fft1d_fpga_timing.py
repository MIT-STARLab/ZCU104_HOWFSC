import matplotlib.pyplot as plt
import numpy as np

fft_size = [128, 256, 512, 1024, 2048]
fpga_hw = [0.000039928, 0.000039756, 0.000055582, 0.000069563, 0.000097487]
fpga_latency = [559 * 10e-9, 949 * 10e-9, 1728 * 10e-9, 3270 * 10e-9, 6353 * 10e-9]
fpga_cosim_time = [868 * 10e-9, 1381 * 10e-9, 2416 * 10e-9, 4477 * 10e-9, 8581 * 10e-9]
zcu_cpu = [0.000156628, 0.000354608, 0.000785367, 0.0017288, 0.00378292]
mac = [0.0000382, 0.0000866, 0.0001939, 0.000587, 0.0009699]

# Convert time from seconds to milliseconds
fpga_hw_ms = np.array(fpga_hw) * 1000
fpga_latency_ms = np.array(fpga_latency) * 1000
fpga_cosim_ms = np.array(fpga_cosim_time) * 1000
zcu_cpu_ms = np.array(zcu_cpu) * 1000
mac_ms = np.array(mac) * 1000

# First Plot: No Log Scaling with FPGA Hardware, Latency, Co-Simulation Time, and Hardware + Co-Simulation Time
plt.figure(figsize=(12, 8))
plt.plot(fft_size, fpga_hw_ms, 'o-', label="ZCU104 FPGA Hardware Run Time", color='blue', linewidth=2)
plt.plot(fft_size, fpga_latency_ms, 'x--', label="ZCU104 FPGA Kernel Latency", color='cyan', linewidth=2)
plt.plot(fft_size, fpga_cosim_ms, 's--', label="ZCU104 FPGA Kernel Co-Simulation Time", color='magenta', linewidth=2)
plt.plot(fft_size, fpga_hw_ms - fpga_cosim_ms, 'd--', label="FPGA Hardware - Co-Simulation Time", color='orange', linewidth=2)

plt.xlabel("FFT Size")
plt.ylabel("Time (milliseconds)")
plt.title("1D FFT Timing for FPGA Implementations")
plt.legend()
plt.grid(True, linestyle='--', linewidth=0.5)
plt.savefig('fft_timing_plot_fpga.png', format='png', dpi=300)
plt.show()

# Second Plot: No Log Scaling with FPGA Hardware, MAC, and CPU
plt.figure(figsize=(12, 8))
plt.plot(fft_size, fpga_hw_ms, 'o-', label="ZCU104 FPGA Hardware Run Time", color='blue', linewidth=2)
plt.plot(fft_size, mac_ms, '^-', label="MAC - Tukey algorithm C++", color='green', linewidth=2)
plt.plot(fft_size, zcu_cpu_ms, 's-', label="ZCU104 CPU - Tukey algorithm C++", color='red', linewidth=2)

plt.xscale('log', base=2)  # Log scale for FFT size
plt.yscale('log')          # Log scale for time

plt.xlabel("Log2(FFT Size)")
plt.ylabel("Time (log scaled) [ms]")
plt.title("1D FFT Timing for FPGA and CPU Implementations (Non-Log Scale)")
plt.legend()
plt.grid(True, linestyle='--', linewidth=0.5)
plt.savefig('fft_timing_plot_fpga_mac_cpu.png', format='png', dpi=300)
plt.show()
