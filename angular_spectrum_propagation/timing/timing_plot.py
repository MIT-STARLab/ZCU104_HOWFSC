import matplotlib.pyplot as plt
import numpy as np

# Data
arr_size = [128, 256, 512, 1024, 2048]

fpga_hw_arth = [0,0,0,0.023984,0]
fpga_hw_geom = [0,0,0,0.023984,0]

zcu_cortex_a53_cpp_arth = [0,0,0,11.302552,0]
zcu_cortex_a53_cpp_geom = [0,0,0,11.302552,0]

mac_m1_cpp_arth = [0.0136593, 0.0632588, 0.273278, 1.1969, 7.53648]
mac_m1_cpp_geom = [0.0135389, 0.062242, 0.272312, 1.19662, 7.53569]

mac_m1_py_arth = [0.001685, 0.004513, 0.016742, 0.071771, 0.305984]
mac_m1_py_geom = [0.001234, 0.004483, 0.016649, 0.071624, 0.305663]

# Convert time from seconds to milliseconds
fpga_hw_arth = np.array(fpga_hw_arth) * 1000
fpga_hw_geom = np.array(fpga_hw_geom) * 1000
zcu_cortex_a53_cpp_arth = np.array(zcu_cortex_a53_cpp_arth) * 1000
zcu_cortex_a53_cpp_geom = np.array(zcu_cortex_a53_cpp_geom) * 1000
mac_m1_cpp_arth = np.array(mac_m1_cpp_arth) * 1000
mac_m1_cpp_geom = np.array(mac_m1_cpp_geom) * 1000
mac_m1_py_arth = np.array(mac_m1_py_arth) * 1000
mac_m1_py_geom = np.array(mac_m1_py_geom) * 1000





# Plot
plt.figure(figsize=(12, 8))  # Increase figure size for better readability

# Plot with improved styles
plt.plot(arr_size, mac_m1_cpp_geom, '^-', label="MAC M1 - C++", color='green', markersize=10, linewidth=2)
plt.plot(arr_size, mac_m1_py_geom, 'o-', label="MAC M1 - Numpy",color='blue', markersize=10, linewidth=2)
plt.plot(arr_size, fpga_hw_geom, 's--', label="ZCU104 FPGA Hardware Run Time", color='orange', markersize=10, linewidth=2)
plt.plot(arr_size, zcu_cortex_a53_cpp_geom, 'x-', label="ZCU104 MPCPU Cortex a53 - C++", color='red', markersize=8, linewidth=2)

# Plot settings
plt.xscale('log', base=2)  # Log scale for matrix size
plt.yscale('log')          # Log scale for time
plt.xlabel("Wavefront Array Size (Log2 scale)", fontsize=14)
plt.ylabel("Time [ms] (Log scale)", fontsize=14)
plt.title("Angular Spectrum Method Timing Analysis for FPGA and CPU Implementations", fontsize=16, fontweight='bold')
plt.legend(loc='best', fontsize=12)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()  # Adjust plot to fit labels and title

# Save and show
# plt.savefig('ang_spec_timing_plot_fpga_mac_cpu.png', format='png', dpi=300)
plt.show()
