import matplotlib.pyplot as plt
import numpy as np

fft_size = [128, 256, 512, 1024, 2048]
fpga = [0.000039928, 0.000039756, 0.000055582, 0.000069563, 0.000097487]
zcu_cpu = [0.000156628, 0.000354608, 0.000785367, 0.0017288, 0.00378292]
mac = [0.0000382, 0.0000866, 0.0001939, 0.000587, 0.0009699]

# Convert time from seconds to milliseconds
fpga_ms = np.array(fpga) * 1000
zcu_cpu_ms = np.array(zcu_cpu) * 1000
mac_ms = np.array(mac) * 1000

plt.figure(figsize=(10, 6))

plt.plot(fft_size, fpga_ms, 'o-', label="ZCU104 FPGA IP Core", color='blue', linewidth=2)
plt.plot(fft_size, zcu_cpu_ms, 's-', label="ZCU104 CPU - Tukey algorithm C++", color='red', linewidth=2)
plt.plot(fft_size, mac_ms, '^-', label="MAC - Tukey algorithm C++", color='green', linewidth=2)

plt.xscale('log', base=2)  # Log scale for FFT size
plt.yscale('log')  # Log scale for time
plt.xlabel("FFT Size")
plt.ylabel("Time (milliseconds)")
plt.title("Timing of 1D FFT for Different Implementations")
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Save the plot as a PNG file
plt.savefig('fft_timing_plot.png', format='png', dpi=300)

# Show the plot
plt.show()
