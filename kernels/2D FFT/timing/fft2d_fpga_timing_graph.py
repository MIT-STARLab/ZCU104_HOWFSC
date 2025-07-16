import matplotlib.pyplot as plt
import numpy as np

# Data
fft_size = np.array([128, 256, 512, 1024, 2048])
fpga_hw = np.array([0.000311, 0.00103, 0.003753, 0.015089, 0.065435])
fft_np = np.array([0.000156, 0.000562, 0.002917, 0.013694, 0.065932])
zcu_cpu = np.array([0.041216, 0.186729, 1.144043, 5.142044, 23.538236])
mac = np.array([0.0066333, 0.0295523, 0.135167, 0.614293, 3.92586])

# Convert time from seconds to milliseconds
fpga_hw_ms = np.array(fpga_hw) * 1000
fft_np_ms = np.array(fft_np) * 1000
zcu_cpu_ms = np.array(zcu_cpu) * 1000
mac_ms = np.array(mac) * 1000


# Convert to log base 2
log_n = np.log2(fft_size)
log_time = np.log(fpga_hw)  # natural log, or use np.log2(fpga_time_sec) if you want base-2

# Linear fit: log(time) = a * log2(n) + b
a, b = np.polyfit(log_n, log_time, 1)

print("2D FFT: log(time) = a * log2(n) + b")
print(f"a (slope): {a:.3f}")
print(f"b (intercept): {b:.3f}")


print(fft_np/fpga_hw)
print(mac/fpga_hw)
print(zcu_cpu/fpga_hw)



# # Plot
# plt.figure(figsize=(12, 8))  # Increase figure size for better readability

# # Plot with improved styles

# plt.plot(fft_size, mac_ms, '^-', label="MAC M1 - Tukey Algorithm C++", color='green', markersize=10, linewidth=2)
# plt.plot(fft_size, fft_np_ms, 'o-', label="MAC M1 - Numpy",color='blue', markersize=10, linewidth=2)
# plt.plot(fft_size, fpga_hw_ms, 's--', label="ZCU104 FPGA Hardware Run Time", color='orange', markersize=10, linewidth=2)
# plt.plot(fft_size, zcu_cpu_ms, 'x-', label="ZCU104 MPCPU Cortex a53 - Tukey Algorithm C++", color='red', markersize=8, linewidth=2)

# # Plot settings
# # plt.xscale('log', base=2)  # Log scale for FFT size
# # plt.yscale('log')          # Log scale for time
# plt.xlabel("FFT Size", fontsize=20)
# plt.ylabel("Time [ms]", fontsize=20)
# plt.title("2D FFT Timing Analysis", fontsize=24, fontweight='bold')

# plt.tick_params(axis='both', which='major', labelsize=16, width=2, length=6, color='black', labelcolor='black')
# plt.tick_params(axis='both', which='minor', labelsize=16, width=2, length=4, color='black', labelcolor='black')

# plt.legend(loc='best', fontsize=20)
# plt.grid(True, which='both', linestyle='--', linewidth=0.5)
# plt.tight_layout()  # Adjust plot to fit labels and title

# # Save and show
# plt.savefig('fft2d_timing_plot_fpga_mac_cpu_real.png', format='png', dpi=300)
# # plt.show()
